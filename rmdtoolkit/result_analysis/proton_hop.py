# Python 3.6.1

import os

import numpy as np

from rmdtoolkit.basic.result import Result


class ProtonHop(Result):
    def __init__(self):
        super().__init__()
        self.delta_t = 1
        self.h_list = list()
        self.t_list = list()
        self.pivot_list = list()

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'DELTA_T':  # unit of this is 'frame' (each frame occupy 1 timestep in LAMMPS)
                            self.delta_t = int(line[1])
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def calculate_forward_hop(self):
        os.chdir(self.work_dir)
        last_frame = 0
        self.find_file(trj_flag=False)
        self.open_file(trj_flag=False)
        while True:
            self.evb_read_frame()
            if not self.evb_info:
                break
            if self.evb_time == 'EOF':
                break
            if self.evb_time < last_frame:
                continue
            last_frame = self.evb_time

            self.pivot_list.append(np.array(self.evb_info['ReactionCenters'])[:, 1])
        self.close_file(trj_flag=False)

        # donor_ids = [0] * len(pivot_ids_list[0])
        # for i, pivot_ids in enumerate(pivot_ids_list):
        #     for j, pivot_id in enumerate(pivot_ids):
        #         if i == 0:
        #             hs = np.zeros(len(pivot_ids_list[0]))
        #         elif i == 1:
        #             hs[j] += 1
        #             donor_ids[j] = pivot_ids_list[i-1][j]
        #         elif pivot_id != pivot_ids_list[i-1][j]:
        #             hs[j] += 1
        #             if pivot_id == donor_ids[j]:
        #                 hs[j] -= 2
        #             donor_ids[j] = pivot_ids_list[i-1][j]
        #         else:
        #             donor_ids[j] = pivot_ids_list[i-1][j]
        #     self.h_list.append(list(hs))  # It is a must to do a conversion. Shallow copy thingy.

    def save_forward_hop_result(self):
        with open('%s%s.hop' % (self.save_dir, self.save_tag), 'w') as save_file:
            for i, (t, hs) in enumerate(zip(self.t_list, self.h_list)):
                save_file.write('{t} {hs} {pivots}\n'.format(t=t, hs=' '.join(list(map(str, hs))),
                                                             pivots=' '.join(list(map(str, self.pivot_list[i])))))
