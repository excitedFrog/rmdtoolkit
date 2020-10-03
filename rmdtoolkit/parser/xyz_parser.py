# Python 3.6.1

import os
import copy
import json
import numpy as np

from .base import ParserBase

THIS_DIR = os.path.dirname(__file__)


class XYZParser(ParserBase):
    def __init__(self, file, process_comment=False, **kwargs):
        super().__init__(file, **kwargs)

        with open('{}/templates/xyz/xyz_frame.json'.format(THIS_DIR), 'r') as fframe:
            self.frame_template = json.load(fframe)

        self.frame_time = 0
        self.frame_json = copy.deepcopy(self.frame_template)
        self.frame_json_list = list()
        self.frame_comment_list = list()
        self.frame_time_list = list()

        self._process_comment = process_comment

    def read_frame(self):
        frame = list()
        status = 0
        at_frame_init = True
        num_atoms = 0

        while True:
            if self.read_line() == 1:
                status = 1  # EOF
                break

            frame.append(self.line)
            if at_frame_init:
                at_frame_init = False
                num_atoms = int(self.line)
            else:
                if len(frame) == num_atoms + 2:
                    break
        self.frame = self.split_list(frame)
        return status

    def parse_frame(self, **kwargs):
        self.frame_json = copy.deepcopy(self.frame_template)

        for i, linel in enumerate(self.frame):
            if not linel:
                continue

            if i == 0:
                self.frame_json["num_atoms"] = int(linel[0])
            elif i == 1:
                comment = ' '.join(linel)
                self.frame_json["comment"] = ' '.join(linel)
                self.frame_comment_list.append(comment)
            else:
                info = np.array(self.frame[i:])
                self.frame_json["atoms"] = info[:, 0].astype(str)
                self.frame_json["coords"] = info[:, 1:4].astype(float)
                break

    def parse_file(self, **kwargs):
        while True:
            status = self.read_frame()
            if status == 1:
                break
            self.frame_time += 1
            self.time_tell()
            self.parse_frame(**kwargs)
            self.frame_json_list.append(copy.deepcopy(self.frame_json))

