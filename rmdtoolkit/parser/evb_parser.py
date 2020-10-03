# Python 3.6.1

import os
import copy
import json
import numpy as np

from .base import ParserBase

THIS_DIR = os.path.dirname(__file__)


class EVBParser(ParserBase):
    def __init__(self, file, **kwargs):
        super().__init__(file, **kwargs)

        with open('{}/templates/evb/evb_frame.json'.format(THIS_DIR), 'r') as fframe, \
                open('{}/templates/evb/evb_complex.json'.format(THIS_DIR), 'r') as fcomplex, \
                open('{}/templates/evb/evb_complex_state.json'.format(THIS_DIR), 'r') as fstate:
            self.frame_template = json.load(fframe)
            self.complex_template = json.load(fcomplex)
            self.state_template = json.load(fstate)

        self.frame_time = int()
        self.frame_json = copy.deepcopy(self.frame_template)
        self.frame_time_list = list()
        self.frame_json_list = list()

    def read_frame(self):
        frame = list()
        status = 0
        break_trigger = 0  # +1 every time reading '***', break at 2

        while True:
            if self.read_line() == 1:
                status = 1
                break
            if self.line.startswith('*****'):
                break_trigger += 1
            if break_trigger > 1:
                self.file.seek(self.last_pos)
                break
            frame.append(self.line)
        self.frame = self.split_list(frame)
        return status

    def parse_frame(self, complex_type=False):
        self.frame_json = copy.deepcopy(self.frame_template)
        icomp = -1

        for i, linel in enumerate(self.frame):
            if not linel:
                continue

            if linel[0] == 'TIMESTEP':
                self.frame_json["time"] = int(linel[1])
            if linel[0] == 'COMPLEX_COUNT':
                self.frame_json["complex_count"] = int(linel[1])
                for j in range(int(linel[1])):
                    self.frame_json["complex"].append(copy.deepcopy(self.complex_template))
            if linel[0] == 'REACTION_CENTER_LOCATION':
                for j in range(i + 1, i + 1 + self.frame_json["complex_count"]):
                    self.frame_json["rc_location"].append(int(self.frame[j][1]))
            if linel[0] == 'ENE_ENVIRONMENT':
                self.frame_json["energy_summary"]["ene_environment"]["total"] = float(linel[1])
            if linel[0].startswith('ENE_COMPLEX'):
                self.frame_json["energy_summary"]["ene_complex"].append(float(linel[1]))
            if linel[0] == 'ENE_INTER_CPLX':
                self.frame_json["energy_summary"]["ene_inter_cplx"] = float(linel[1])
            if linel[0] == 'ENE_TOTAL':
                self.frame_json["energy_summary"]["ene_total"] = float(linel[1])
            if linel[0] == 'ENVIRONMENT':
                keys = ''.join(linel[1:]).replace('[', '').replace(']', '').split('|')
                for key, value in zip(keys, self.frame[i+1]):
                    self.frame_json["energy_summary"]["ene_environment"][key] = float(value)

            if linel[0] == 'START_OF_COMPLEX':
                icomp += 1
            if linel[0] == 'COMPLEX' and linel[1] == str(icomp+1) + ':':
                self.frame_json["complex"][icomp]["states_count"] = int(linel[2])
                for j in range(int(linel[2])):
                    self.frame_json["complex"][icomp]["states"].append(copy.deepcopy(self.state_template))
            if linel[0] == 'STATES':
                keys = ''.join(linel[1:]).replace('[', '').replace(']', '').split('|')
                for istate, j in enumerate(range(i + 1, i + 1 + self.frame_json["complex"][icomp]["states_count"])):
                    for key, value in zip(keys, self.frame[j]):
                        self.frame_json["complex"][icomp]["states"][istate][key] = int(value)
            if linel[0] == 'EIGEN_VECTOR':
                self.frame_json["complex"][icomp]["eigen_vector"] = list(map(float, self.frame[i + 1]))
            if linel[0] == 'CEC_COORDINATE':
                self.frame_json["complex"][icomp]["cec_coordinate"] = list(map(float, self.frame[i + 1]))

        self.frame_time = self.frame_json["time"]

        if complex_type:  # "1-${NUM_SHELL_1}-${NUM_SHELL_2}-..."
            self.frame_json["extra"]["complex_type"] = list()
            for comp in self.frame_json["complex"]:
                shells = np.array([])
                mol_bs = np.array([])
                num_shells = list()
                for state in comp["states"]:
                    shells = np.append(shells, state["shell"])
                    mol_bs = np.append(mol_bs, state["mol_B"])
                for i_shell in np.unique(shells):
                    num_shell = len(mol_bs[shells == i_shell])
                    num_shells.append(str(num_shell))
                comp_type = '-'.join(num_shells)
                self.frame_json["extra"]["complex_type"].append(comp_type)

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
