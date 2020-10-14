from datetime import timedelta
from glob import glob
import os
from os import PathLike
from pathlib import Path
import re
import shutil
import tarfile

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.server import SlurmConfig
from nemspy import ModelingSystem
from nemspy.model import ADCIRCEntry
import numpy
import requests

from .job_script import EnsembleSlurmScript, SlurmEmailType
from ..utilities import get_logger, repository_root

LOGGER = get_logger('configuration.adcirc')

TACC_TASKS_PER_NODE = 68


def download_test_configuration(directory: str):
    """
    fetch shinnecock inlet test data
    :param directory: local directory
    """

    if not isinstance(directory, Path):
        directory = Path(directory)

    if not directory.exists():
        os.makedirs(directory, exist_ok=True)

    url = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'
    remote_file = requests.get(url, stream=True)
    temporary_filename = directory / 'temp.tar.gz'
    with open(temporary_filename, 'b+w') as local_file:
        local_file.write(remote_file.raw.read())
    with tarfile.open(temporary_filename, 'r:bz2') as local_file:
        local_file.extractall(directory)
    os.remove(temporary_filename)


def write_adcirc_configurations(
        nems: ModelingSystem,
        runs: {str: (float, str)},
        input_directory: PathLike,
        output_directory: PathLike,
        name: str = None,
        partition: str = None,
        email_address: str = None,
        tacc: bool = False,
        wall_clock_time: timedelta = None,
):
    """
    Generate ADCIRC run configuration for given variable values.

    :param runs: dictionary of run name to run value and mesh attribute name
    :param nems: NEMSpy ModelingSystem object, populated with models and connections
    :param input_directory: path to input data
    :param output_directory: path to store run configuration
    :param name: name of this perturbation
    :param partition: Slurm partition
    :param email_address: email address
    :param tacc: whether to configure for TACC
    :param wall_clock_time: wall clock time of job
    """

    if not isinstance(input_directory, Path):
        input_directory = Path(input_directory)
    if not isinstance(output_directory, Path):
        input_directory = Path(output_directory)

    if not input_directory.exists():
        os.makedirs(input_directory, exist_ok=True)
    if not output_directory.exists():
        os.makedirs(output_directory, exist_ok=True)

    if name is None:
        name = 'perturbation'

    if 'ocn' not in nems or not isinstance(nems['ocn'], ADCIRCEntry):
        nems['ocn'] = ADCIRCEntry(11)

    fort14_filename = input_directory / 'fort.14'
    nems_executable = output_directory / 'NEMS.x'

    launcher = 'ibrun' if tacc else 'srun'
    run_name = 'ADCIRC_GAHM_GENERIC'

    if partition is None:
        partition = 'development'

    if wall_clock_time is None:
        wall_clock_time = timedelta(minutes=30)

    if not fort14_filename.is_file():
        download_test_configuration(input_directory)

    # open mesh file
    mesh = AdcircMesh.open(fort14_filename, crs=4326)

    # init tidal forcing and setup requests
    tidal_forcing = Tides()
    tidal_forcing.use_all()

    mesh.add_forcing(tidal_forcing)

    slurm = SlurmConfig(
        account=None,
        ntasks=nems.processors,
        run_name=run_name,
        partition=partition,
        walltime=wall_clock_time,
        nodes=divmod(nems.processors, TACC_TASKS_PER_NODE)[0] if tacc else None,
        mail_type='all' if email_address is not None else None,
        mail_user=email_address,
        log_filename=f'{name}.log',
        modules=['intel', 'impi', 'netcdf'],
        path_prefix='$HOME/adcirc/build',
        launcher=launcher,
    )

    # instantiate AdcircRun object.
    driver = AdcircRun(
        mesh=mesh,
        start_date=nems.start_time,
        end_date=nems.start_time + nems.duration,
        spinup_time=timedelta(days=5),
        server_config=slurm,
    )
    driver.import_stations(Path(repository_root()) / 'examples/data/stations.txt')
    driver.set_elevation_stations_output(timedelta(minutes=6), spinup=timedelta(minutes=6))
    driver.set_elevation_surface_output(timedelta(minutes=6), spinup=timedelta(minutes=6))
    driver.set_velocity_stations_output(timedelta(minutes=6), spinup=timedelta(minutes=6))
    driver.set_velocity_surface_output(timedelta(minutes=6), spinup=timedelta(minutes=6))

    for run_name, (value, attribute_name) in runs.items():
        run_directory = output_directory / run_name
        LOGGER.info(f'writing config to "{run_directory}"')
        if not isinstance(value, numpy.ndarray):
            value = numpy.full([len(driver.mesh.coords)], fill_value=value)
        if not driver.mesh.has_attribute(attribute_name):
            driver.mesh.add_attribute(attribute_name)
        driver.mesh.set_attribute(attribute_name, value)
        driver.write(run_directory, overwrite=True)
        nems.write(run_directory, overwrite=True)
        if nems_executable.exists():
            shutil.copyfile(nems_executable, run_directory / 'NEMS.x')

    ensemble_slurm_script = EnsembleSlurmScript(
        account=None,
        tasks=nems.processors,
        duration=wall_clock_time,
        partition=partition,
        launcher=launcher,
        run=run_name,
        email_type=SlurmEmailType.ALL if email_address is not None else None,
        email_address=email_address,
        log_filename=f'{name}.log',
        nodes=divmod(nems.processors, TACC_TASKS_PER_NODE)[0] if tacc else None,
        modules=['intel', 'impi', 'netcdf'],
        path_prefix='$HOME/adcirc/build',
    )
    ensemble_slurm_script.write(output_directory, overwrite=True)

    pattern = re.compile(' p*adcirc')
    replacement = ' NEMS.x'
    for job_filename in glob(str(output_directory / '**' / 'slurm.job'), recursive=True):
        with open(job_filename) as job_file:
            text = job_file.read()
        matched = pattern.search(text)
        if matched:
            LOGGER.debug(
                f'replacing `{matched.group(0)}` with `{replacement}`' f' in "{job_filename}"'
            )
            text = re.sub(pattern, replacement, text)
            with open(job_filename, 'w') as job_file:
                job_file.write(text)
