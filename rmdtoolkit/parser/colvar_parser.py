# Python 3.6.1

import os
import json
import numpy as np

from .base import ParserBase

THIS_DIR = os.path.dirname(__file__)


class COLVARParser(ParserBase):
    def __init__(self, file, time_col=0, **kwargs):
        super().__init__(file, **kwargs)

        with open('{}/templates/colvar/colvar_file.json'.format(THIS_DIR), 'r') as template:
            self.file_json = json.load(template)

        self.frame_time_list = list()

        self._time_col_i = time_col

    def parse_file(self, start=0, stop=np.inf, t_mult=1):
        array = list()
        for i, line in enumerate(self.file):
            if i == 0:
                ll = line.split()
                self.file_json['col_names'] = ll[2:]
            else:
                if line.startswith('#!'):
                    continue
                ll = list(map(float, line.strip().split()))
                time = ll[self._time_col_i] * t_mult
                if time < start:
                    continue
                if time > stop:
                    break
                array.append(ll)
                self.frame_time_list.append(time)
        array = np.array(array).T
        self.file_json['col_vectors'] = array

    def read_frame(self, **kwargs):
        raise NotImplementedError

    def parse_frame(self, **kwargs):
        raise NotImplementedError
