# Python 3.6.1

import sys

from ..config import ConfigManager


class InputBase(object):
    def __init__(self, work_dir, base_name,
                 err=sys.stderr, out=sys.stdout,
                 config_path=None):
        self.cm = ConfigManager(config_path=config_path)

        self.work_dir = work_dir
        self.base_name = base_name

        self.err = err
        self.out = out
