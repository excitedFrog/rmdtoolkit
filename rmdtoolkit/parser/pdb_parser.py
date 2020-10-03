# Python 3.6.1

import os
import copy
import json
import numpy as np

from .base import ParserBase

THIS_DIR = os.path.dirname(__file__)


class PDBParser(ParserBase):
    def __init__(self, file, **kwargs):
        super().__init__(file, **kwargs)
        """
        :param path: absolute path of the file to be parsed
        :param as_df: return as Pandas DataFrame or not
        """

        with open('{}/templates/pdb/pdb_frame.json'.format(THIS_DIR), 'r') as fframe:
            self.frame_template = json.load(fframe)

        self.frame_json = copy.deepcopy(self.frame_template)

    def read_frame(self, **kwargs):
        self.read_lines()

    def parse_frame(self, **kwargs):
        self.frame = list(map(lambda string: string.strip('\n'), self.frame))
        self.frame = np.array(self.frame)[list(map(lambda string: string.startswith('ATOM'), self.frame))]

        self.frame_json["atoms"] = np.array(list(map(lambda string: string[13:16].strip(' '), self.frame)))
        self.frame_json["elements"] = np.array(list(map(lambda string: string[76:78].strip(' '), self.frame)))
        self.frame_json["res_names"] = np.array(list(map(lambda string: string[17:20].strip(' '), self.frame)))
        self.frame_json["chain_ids"] = np.array(list(map(lambda string: string[21:22].strip(' '), self.frame)))
        self.frame_json["res_seqs"] = np.array(list(map(lambda string: string[22:26].strip(' '), self.frame)))
        self.frame_json["coords"] = np.array(list(map(lambda string: list(map(float, string[30:54].split())),
                                                      self.frame)))

    def parse_file(self, **kwargs):
        self.read_frame()
        self.parse_frame()
