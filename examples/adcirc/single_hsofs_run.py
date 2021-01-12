#! /usr/bin/env python

from datetime import datetime, timedelta
from pathlib import Path
import sys

from nemspy import ModelingSystem
from nemspy.model import ADCIRCEntry, AtmosphericMeshEntry, WaveMeshEntry

sys.path.append(str(Path(__file__).absolute().parent.parent.parent))

from ensemble_perturbation.configuration.adcirc import download_shinnecock_mesh, write_adcirc_configurations
from ensemble_perturbation.configuration.job_script import HPC
from ensemble_perturbation.utilities import repository_root

DATA_DIRECTORY = repository_root() / 'examples/data'
#INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'hsofs'
INPUT_DIRECTORY = Path("/scratch2/COASTAL/coastal/save/Saeed.Moghimi/setups/nems_inp/hsofs_grid_v1/")
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'configuration' / 'hsofs'

if __name__ == '__main__':
    runs = {f'nems_hsofs_test': (None, None)}

    if not (INPUT_DIRECTORY / 'fort.14').exists():
        download_shinnecock_mesh(INPUT_DIRECTORY)

    nems = ModelingSystem(
        start_time=datetime(2012, 10, 22, 6),
        duration=timedelta(days=14.5),
        interval=timedelta(hours=1),
        atm=AtmosphericMeshEntry('/scratch2/COASTAL/coastal/save/Saeed.Moghimi/setups/nems_inp/hsofs_forcings/san_v2/inp_atmesh/SANDY_HWRF_HSOFS_Nov2018.nc'),
        wav=WaveMeshEntry('/scratch2/COASTAL/coastal/save/Saeed.Moghimi/setups/nems_inp/hsofs_forcings/san_v2/inp_ww3/ww3.HWRF.NOV2018.2012_sxy.nc'),
        ocn=ADCIRCEntry(382),
    )

    nems.connect('ATM', 'OCN')
    nems.connect('WAV', 'OCN')
    nems.sequence = [
        'ATM -> OCN',
        'WAV -> OCN',
        'ATM',
        'WAV',
        'OCN',
    ]

    write_adcirc_configurations(
        nems,
        runs,
        INPUT_DIRECTORY,
        OUTPUT_DIRECTORY,
        name='nems_hsofs_test',
        email_address='zachary.burnett@noaa.gov',
        platform=HPC.HERA,
        spinup=timedelta(days=12.5),
    )
