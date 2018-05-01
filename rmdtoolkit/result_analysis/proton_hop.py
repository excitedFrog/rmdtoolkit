# Python 3.6.1

import os

import numpy as np
from itertools import accumulate as acc

from rmdtoolkit.basic.result import Result


class ProtonHop(Result):
    def __init__(self):
        super().__init__()
        self.h_lists = list()
        self.pivot_lists = list()  # Also good for multiple CECs
        self.time_list = list()

    def read_input(self):
        super().read_input()

    def forward_hop_worker(self):
        self.trj_flag = False
        self.evb_flag = True
        self.input_flag = False
        self.analysis_template(inner_compute_func=self.log_forward_hop_info,
                               inner_save_func=self.void_func,
                               outer_compute_func=self.calculate_forward_hop,
                               outer_save_func=self.save_forward_hop_result)

    def log_forward_hop_info(self):
        self.pivot_lists.append(np.array(self.evb_info['ReactionCenters'])[:, 1])
        self.time_list.append(self.evb_time)

    def calculate_forward_hop(self):
        self.pivot_lists = np.array(self.pivot_lists).T
        # I really hate this chunk of for-loop and if/elif/else, but for the sake of readability I'll keep it this way.
        rattle_trigger = False
        for pivot_list in self.pivot_lists:
            h_list = list()
            for i, pivot_id in enumerate(pivot_list):
                if i == 0:
                    delta_h = 0
                elif not pivot_id - pivot_list[i-1]:
                    delta_h = 0
                elif pivot_id - pivot_list[i-2]:
                    delta_h = 1
                    rattle_trigger = False
                elif rattle_trigger:
                    delta_h = 1
                    rattle_trigger = False
                else:
                    delta_h = -1
                    rattle_trigger = True
                h_list.append(delta_h)
            h_list = list(acc(h_list))
            self.h_lists.append(h_list)

    def save_forward_hop_result(self):
        self.time_list = np.array(self.time_list).astype(str)
        self.h_lists = np.array(self.h_lists).astype(str)
        result_arr = np.vstack((self.time_list, self.h_lists)).T
        with open('{}{}.hop'.format(self.save_dir, self.save_tag), 'w') as save_file:
            for result_line in result_arr:
                save_file.write(' '.join(result_line) + '\n')
