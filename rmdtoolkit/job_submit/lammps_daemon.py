# Python 3.6.1

import subprocess
import os
import time
import sys
import signal

from job_submit.daemon import Daemon
from config.config_manager import ConfigManager

CM = ConfigManager()
TMP_DIR = CM.get('job_daemon', 'daemon_tmp_dir')
CURRENT_DIR = os.getcwd()


class LDaemon(Daemon):
    def __init__(self, daemon_pidfile='%sdaemon.pid/' % TMP_DIR, child_pidfile='%schild.pid/' % TMP_DIR,
                 stdin='/dev/null', stdout='/dev/null', stderr='/dev/null'):
        super().__init__(daemon_pidfile=daemon_pidfile, stdin=stdin, stdout=stdout, stderr=stderr)
        self.child_pidfile = child_pidfile
        self.process = None

    def log_child_pid(self):
        with open(self.child_pidfile, 'w') as f:
            print(self.process.pid, file=f)

    def run(self, restart_flag=False):
        sys.stdout.write('Daemon started with pid {}\n'.format(os.getpid()))
        sys.stdout.flush()
        if restart_flag:
            self.process = subprocess.Popen(['python', 'job_monitor_run.py', '-r'])
        else:
            self.process = subprocess.Popen(['python', 'job_monitor_run.py'])
        sys.stdout.write('Child spawned! {}\n'.format(self.process.pid))
        sys.stdout.flush()
        self.log_child_pid()
        while True:
            time.sleep(60)
            sys.stdout.write('Daemon Alive! {}\n'.format(time.ctime()))
            sys.stdout.flush()
            if self.process.poll() is not None:
                sys.stdout.write('Child died!\n')
                sys.stdout.flush()
                self.process = subprocess.Popen(['python', '{}/zLAMMPS.py'.format(CURRENT_DIR), '-r'])
                sys.stdout.write('Child respawned! {}\n'.format(self.process.pid))
                sys.stdout.flush()
                self.log_child_pid()

    def start(self, restart_flag=False):
        try:
            self.daemonize()
        except RuntimeError as error:
            print(error, file=sys.stderr)
            raise SystemExit(1)
        self.run(restart_flag=restart_flag)

    def stop(self):
        super().stop()
        with open(self.child_pidfile, 'r') as f:
            child_pid = int(f.read())
        os.kill(child_pid, signal.SIGKILL)
        os.remove(self.child_pidfile)
