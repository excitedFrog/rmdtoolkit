# Python 3.6.1

import sys

from job_submit.job_monitor import JobMonitor
from config.config_manager import ConfigManager

is_restart = False
if '-r' in sys.argv:
    is_restart = True

CM = ConfigManager()
PREFIX_DIR = CM.get('lammps', 'simulation_prefix_dir')

JM = JobMonitor()
JM.prefix_dir = PREFIX_DIR
for key, item in CM.get_section('iter_dirs'):
    JM.iter_lists.append(item.split())

JM.run(is_restart)
