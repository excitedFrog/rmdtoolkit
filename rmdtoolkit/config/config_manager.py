# Python 3.6.1

import warnings
import configparser
from os import path

here = path.abspath(path.dirname(__file__))


class ConfigManager(object):
    def __init__(self):
        self.config_path = self.get_config_path()
        self.cp = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        self.cp.read(self.config_path)

    def get(self, section, option):
        info = self.cp.get(section, option)
        if info:
            return info
        else:
            raise Exception('[CONFIG ERROR] <{}> in section <{}> is not set!'.format(option, section))

    def set(self, section, option, value):
        if not self.cp.has_section(section):
            warnings.warn('[CONFIG WARNING] section <{}> does not exist. Make sure there is not a typo.'
                          .format(section))
            self.cp.add_section(section)
        if not self.cp.has_option(section, option):
            warnings.warn('[CONFIG WARNING] section <{}> does not have option <{}>. Make sure there is not a typo.'
                          .format(section, option))
        self.cp.set(section, option, value)
        self.cp.write(open(self.config_path, 'w'))
        print('[CONFIG] option <{}> in section <{}> has been set as \'{}\''.format(option, section, value))

    def get_list(self, section, option):
        return self.cp.get(section, option).split()

    def get_section(self, section):
        return self.cp.items(section)

    def print_section(self, section):
        info = self.get_section(section)
        print('[{}]'.format(section))
        print('\n'.join([' = '.join(_) for _ in info]))

    @staticmethod
    def set_config_path(string):
        with open(path.join(here, 'config_path'), 'w', encoding='utf-8') as config_path_file:
            config_path_file.write(string)

    @staticmethod
    def get_config_path():
        with open(path.join(here, 'config_path'), 'r', encoding='utf-8') as config_path_file:
            config_path = config_path_file.readline()
        return config_path
