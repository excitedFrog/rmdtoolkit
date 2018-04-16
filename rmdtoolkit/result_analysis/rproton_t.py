# Python 3.6.1

import os
from collections import OrderedDict

import numpy as np

from basic.result import Result


class RprotonT(Result):
    def __init__(self):
        super().__init__()
        self.rt_input_file = None
        self.stop_frame = float('inf')
        self.frame_tot = 0

        self.specie_tags = list()
        self.specie_numbers = list()
        self.total_species = 0
        self.sections = OrderedDict()

        self.sys_names = list()
        self.r_cec_list = list()

    def read_input(self):
        super().read_input()
        with open(self.rt_input_file, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'STOP_FRAME':  # stop frame (unit is frame)
                            self.stop_frame = int(line[1])
                        elif line[0] == 'SYSNAMES':  # should be in order
                            self.sys_names = line[1:]
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def calculate_rt_cec(self):
        os.chdir(self.work_dir)
        if not self.sys_names:
            self.sys_names = sorted([file_name.strip('.lammpstrj')
                                     for file_name in os.listdir(self.work_dir) if file_name.endswith('.lammpstrj')])
        is_first_frame = True
        system_size = 0
        last_frame = 0
        for sys_name in self.sys_names:
            self.sys_name = sys_name
            self.find_file()
            self.open_file()
            while True:
                self.trj_read_frame()
                if not is_first_frame and len(self.trj_info) != system_size:
                    break
                if self.trj_time == 'EOF':
                    break
                if self.frame_tot > self.stop_frame:
                    break
                elif self.trj_time < last_frame:
                    continue
                self.evb_match_trj()
                if not self.evb_info:
                    break
                last_frame = self.trj_time
                self.frame_tot += 1

                if is_first_frame:
                    is_first_frame = False
                    system_size = len(self.trj_info)

                if self.frame_tot % 10 == 0:
                    print('{} PROCESSING FRAME {}'.format(self.save_tag, self.trj_time))
                com = self.find_com()
                cecs = self.find_cec()[0]
                r_cecs = list()
                for cec in cecs:
                    r_cec = com - cec
                    r_cec = np.linalg.norm(r_cec)
                    r_cecs.append(str(r_cec))
                self.r_cec_list.append(r_cecs)
            self.close_file()

    def save_rt_result(self):
        with open('{}{}.rt'.format(self.save_dir, self.save_tag), 'w') as save_file:
            for i, cec_list in enumerate(self.r_cec_list):
                save_file.write('{} {}\n'.format(i, ' '.join(cec_list)))
