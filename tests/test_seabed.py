import unittest

import numpy
from pyproj import CRS, Transformer

from ensemble_perturbation.inputs.seabed.ngdc import NGDCSeabedDescriptions


class TestSeabedDescriptions(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.crs = CRS.from_epsg(32618)

        self.surveys = NGDCSeabedDescriptions.all_surveys()[:5]

        bounds = numpy.array([[-77, 39], [-75, 40]])
        transformer = Transformer.from_crs(CRS.from_epsg(4326), self.crs)
        self.bounds = numpy.ravel(numpy.stack(transformer.transform(
            bounds[:, 0], bounds[:, 1]), axis=1))

    def test_seabed_descriptions(self):
        seabed = NGDCSeabedDescriptions(bounds=self.bounds,
                                        surveys=self.surveys, crs=self.crs)

        assert seabed.data.shape[0] > 0
        assert seabed.data.shape[1] == 14
        assert len(seabed.descriptions) > 0

    def test_interpolation(self):
        seabed = NGDCSeabedDescriptions(bounds=self.bounds,
                                        surveys=self.surveys, crs=self.crs)

        seabed.interpolate((-76.5, 39.2), crs=CRS.from_epsg(4326))


if __name__ == '__main__':
    unittest.main()
