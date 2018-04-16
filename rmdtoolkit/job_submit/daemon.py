# Python 3.6.1

import os
import sys
import atexit
import signal

from config.config_manager import ConfigManager

CM = ConfigManager()
ROOT_DIR = CM.get('job_daemon', 'daemon_root_dir')
TMP_DIR = CM.get('job_daemon', 'daemon_tmp_dir')


class Daemon:
    def __init__(self, daemon_pidfile='%sdaemon.pid/' % TMP_DIR, stdin='/dev/null', stdout='/dev/null',
                 stderr='/dev/null'):
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        self.daemon_pidfile = daemon_pidfile

    def daemonize(self):
        if os.path.exists(self.daemon_pidfile):
            print('Daemon already running.')
        # First fork (detach from parent)
        try:
            if os.fork() > 0:
                raise SystemExit(0)
        except OSError as error:
            raise RuntimeError('fork #1 failed: {0} ({1})'.format(error.errno, error.strerror))
        os.chdir(ROOT_DIR)
        os.setsid()
        os.umask(0o22)
        # Second fork (relinquish session leadership)
        try:
            if os.fork() > 0:
                raise SystemExit(0)
        except OSError as error:
            raise RuntimeError('fork #2 failed: {0} ({1})\n'.format(error.errno, error.strerror))
        # Flush I/O buffers
        sys.stdout.flush()
        sys.stderr.flush()
        # Replace file descriptors for stdin, stdout, and stderr
        with open(self.stdin, 'rb', 0) as f:
            os.dup2(f.fileno(), sys.stdin.fileno())
        with open(self.stdout, 'ab', 0) as f:
            os.dup2(f.fileno(), sys.stdout.fileno())
        with open(self.stderr, 'ab', 0) as f:
            os.dup2(f.fileno(), sys.stderr.fileno())
        # Write the PID file
        with open(self.daemon_pidfile, 'w') as f:
            print(os.getpid(), file=f)
        # Arrange to have the PID file removed on exit/signal
        atexit.register(lambda: os.remove(self.daemon_pidfile))
        signal.signal(signal.SIGTERM, self.__sigterm_handler)

    # Signal handler for termination (required)
    @staticmethod
    def __sigterm_handler(signo, frame):
        raise SystemExit(1)

    def start(self):
        try:
            self.daemonize()
        except RuntimeError as error:
            print(error, file=sys.stderr)
            raise SystemExit(1)
        self.run()

    def stop(self):
        try:
            if os.path.exists(self.daemon_pidfile):
                with open(self.daemon_pidfile) as f:
                    os.kill(int(f.read()), signal.SIGTERM)
            else:
                print('Not running.', file=sys.stderr)
                raise SystemExit(1)
        except OSError as error:
            if 'No such process' in str(error) and os.path.exists(self.daemon_pidfile):
                os.remove(self.daemon_pidfile)

    def restart(self):
        self.stop()
        self.start()

    def run(self):
        pass
