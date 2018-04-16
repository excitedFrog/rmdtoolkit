# Python 3.6.1

import sys

from rmdtoolkit.config.config_manager import ConfigManager
from rmdtoolkit.job_submit.lammps_daemon import LDaemon

CM = ConfigManager()
LOG_FILE = CM.get('job_daemon', 'daemon_log_file')

MyDaemon = LDaemon(stdout=LOG_FILE, stderr=LOG_FILE)
if len(sys.argv) != 2:
    print('Usage: {} [start|stop|restart]'.format(sys.argv[0]), file=sys.stderr)
    raise SystemExit(1)

if sys.argv[1] == 'start':
    MyDaemon.start(restart_flag=False)
elif sys.argv[1] == 'stop':
    MyDaemon.stop()
elif sys.argv[1] == 'restart':
    MyDaemon.start(restart_flag=True)
else:
    print('Unknown command {!r}'.format(sys.argv[1]), file=sys.stderr)
    raise SystemExit(1)
