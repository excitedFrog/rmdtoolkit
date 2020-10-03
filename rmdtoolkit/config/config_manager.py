# Python 3.6.1

import os
import warnings
import configparser

THIS_DIR = os.path.dirname(__file__)


class ConfigManager(object):
    def __init__(self, config_path=None):
        if config_path:
            self.set_config_path(config_path)
        self.config_path = self.get_config_path()
        self.cp = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        self.cp.read(self.config_path)

    def get(self, section, option):
        info = self.cp.get(section, option)
        if info:
            return info
        else:
            raise ValueError('[CONFIG ERROR] <{}> in section <{}> is not set!'.format(option, section))

    def set(self, section, option, value):
        if not self.cp.has_section(section):
            warnings.warn('[CONFIG WARNING] section <{}> does not exist.'.format(section))
            self.cp.add_section(section)
        if not self.cp.has_option(section, option):
            warnings.warn('[CONFIG WARNING] section <{}> does not have option <{}>.'.format(section, option))
        self.cp.set(section, option, value)
        self.cp.write(open(self.config_path, 'w'))
        print('[CONFIG] option <{}> in section <{}> has been set as \'{}\''.format(option, section, value), )

    def get_list(self, section, option):
        return self.cp.get(section, option).split()

    def get_section(self, section):
        return self.cp.items(section)

    def print_section(self, section):
        info = self.get_section(section)
        print('[{}]'.format(section))
        print('\n'.join([' = '.join(_) for _ in info]))

    @staticmethod
    def set_config_path(path):
        with open(os.path.join(THIS_DIR, 'ConfigPath.config'), 'w', encoding='utf-8') as config_path_file:
            config_path_file.write(path)

    @staticmethod
    def get_config_path():
        with open(os.path.join(THIS_DIR, 'ConfigPath.config'), 'r', encoding='utf-8') as config_path_file:
            config_path = config_path_file.read()
        return config_path
