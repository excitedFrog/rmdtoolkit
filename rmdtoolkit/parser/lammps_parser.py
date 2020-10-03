# Python 3.6.1

import os
import copy
import json
import numpy as np

from .base import ParserBase

THIS_DIR = os.path.dirname(__file__)


class DumpParser(ParserBase):
    def __init__(self, file, **kwargs):
        super().__init__(file, **kwargs)

        with open('{}/templates/lammps/dump_frame.json'.format(THIS_DIR), 'r') as fframe:
            self.frame_template = json.load(fframe)

        self.frame_time = int()
        self.frame_json = copy.deepcopy(self.frame_template)
        self.frame_time_list = list()
        self.frame_json_list = list()

    def read_frame(self):
        frame = list()
        status = 0
        break_trigger = 0  # +1 every time reading 'ITEM: TIMESTEP', break at 2

        while True:
            if self.read_line() == 1:
                status = 1
                break
            if self.line.startswith('ITEM: TIMESTEP'):
                break_trigger += 1
            if break_trigger > 1:
                self.file.seek(self.last_pos)
                break
            frame.append(self.line)
        self.frame = self.split_list(frame)
        return status

    def parse_frame(self, box_lengths=True):
        self.frame_json = copy.deepcopy(self.frame_template)

        for i, linel in enumerate(self.frame):
            if not linel:
                continue

            if linel[0] == 'ITEM:':
                if linel[1] == 'TIMESTEP':
                    self.frame_json["time"] = int(self.frame[i+1][0])
                if linel[1:4] == ['NUMBER', 'OF', 'ATOMS']:
                    self.frame_json["num_atoms"] = int(self.frame[i+1][0])
                if linel[1:3] == ['BOX', 'BOUNDS']:
                    bounds = np.array([])
                    for j in range(i + 1, i + 4):
                        bounds = np.append(bounds, np.array(self.frame[j]).astype(float))
                    bounds = bounds.reshape((3, 2))
                    self.frame_json["bounds"] = bounds
                    if box_lengths:
                        box_lengths = np.abs(bounds[:, 1] - bounds[:, 0])
                        self.frame_json["box_lengths"] = box_lengths
                if linel[1] == 'ATOMS':
                    atom_keys = linel[2:]
                    self.frame_json["atom_keys"] = atom_keys
                    self.frame_json["atom_info"] = np.array(self.frame[i+1:]).astype(float)
                    break

        self.frame_time = self.frame_json["time"]

    def parse_file(self, debug=None, **kwargs):
        while True:
            status = self.read_frame()
            if status == 1:
                break
            self.parse_frame(**kwargs)
            self.frame_time_list.append(self.frame_time)
            self.frame_json_list.append(copy.deepcopy(self.frame_json))
            self.time_tell()
            if debug and len(self.frame_time_list) >= debug:
                break


# TODO
class DataParser(ParserBase):
    def __init__(self, file, **kwargs):
        super().__init__(file, **kwargs)
