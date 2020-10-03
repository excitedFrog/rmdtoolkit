# Python 3.6.1

import sqlite3 as sql
from io import StringIO

from ..config import ConfigManager


class MyDatabase(object):
    def __init__(self, config_path=None):
        self.cm = ConfigManager(config_path=config_path)
        self.connection = None
        self.cursor = None

        self.db_path = self.cm.get('paths', 'top_database')

    def initialize(self):
        self.connection = sql.connect(self.db_path)
        self.cursor = self.connection.cursor()
