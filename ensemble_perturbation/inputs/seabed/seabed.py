from abc import ABC, abstractmethod
import os

from geopandas import GeoDataFrame
import numpy
from pyproj import CRS, Transformer
from scipy.spatial.ckdtree import cKDTree
from shapely.geometry import MultiPoint

from ensemble_perturbation import get_logger

LOGGER = get_logger('seabed')


class SeabedDescriptions(ABC):
    longitude_field = 'Longitude'
    latitude_field = 'Latitude'
    description_field = 'Description'

    def __init__(self, bounds: (float, float, float, float) = None,
                 surveys: [str] = None, crs: CRS = None):
        self.bounds = bounds
        self.__surveys = surveys
        self.crs = CRS.from_user_input(crs) if crs is not None else None

    @classmethod
    @abstractmethod
    def all_surveys(cls) -> [str]:
        raise NotImplementedError

    @property
    def surveys(self) -> [str]:
        if self.__surveys is None:
            self.__surveys = self.__class__.all_surveys()
        return self.__surveys

    def __getitem__(self, survey: str) -> GeoDataFrame:
        raise NotImplementedError

    @property
    @abstractmethod
    def data(self) -> GeoDataFrame:
        raise NotImplementedError

    @property
    @abstractmethod
    def descriptions(self) -> [str]:
        raise NotImplementedError

    def __iter__(self) -> GeoDataFrame:
        for survey in self.surveys:
            yield self[survey]

    def write(self, filename: str, **kwargs):
        drivers = {
            '.csv': 'CSV',
            '.gpkg': 'GPKG',
            '.json': 'GeoJSON',
            '.shp': 'Esri Shapefile',
            '.gdb': 'OpenFileGDB',
            '.gml': 'GML',
            '.xml': 'GML'
        }

        extension = os.path.splitext(filename)[-1]
        kwargs['driver'] = drivers[extension]

        self.data.to_file(filename, **kwargs)

    def interpolate(self, points: [(float, float)], crs: CRS = None) -> ():
        if not isinstance(points, numpy.ndarray):
            points = numpy.array(points)
        if len(points.shape) < 2:
            points = numpy.expand_dims(points, axis=0)
        if crs is not None and crs != self.crs:
            transformer = Transformer.from_crs(crs, self.crs)
            points = numpy.stack(transformer.transform(points[:, 0],
                                                       points[:, 1]), axis=1)

        data_points = numpy.array(list(self.data['geometry'].apply(
            lambda point: (point.x, point.y))))

        tree = cKDTree(data_points)
        distances, indices = tree.query(points, 3)
        nearest_points = numpy.squeeze(data_points[indices])

        interpolated_values = []
        if MultiPoint(points).within(MultiPoint(nearest_points).convex_hull):
            nearest_values = [self.data.cx[point[0],
                                           point[1]][self.description_field]
                              for point in nearest_points]

            for point in points:
                nearest_values
