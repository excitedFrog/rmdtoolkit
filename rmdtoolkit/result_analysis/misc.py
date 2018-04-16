# Python 3.6.1

import os

import numpy as np

from basic.result import Result
from config.config_manager import ConfigManager

CM = ConfigManager()


class Misc(Result):
    def __init__(self):
        super().__init__()

    def fetch_relative_cec_pos(self, stop_at_t=float('inf')):
        os.chdir(self.work_dir)
        self.find_file()
        self.open_file()

        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF' or self.trj_time > stop_at_t:
                break
            self.evb_match_trj()
            if self.trj_time % 1000 == 0:
                print('Done Frame %s' % self.trj_time)

            com = self.find_com()
            cecs = self.evb_info['CECS']
            with open('%s%s.cecvecs' % (self.save_dir, self.save_tag), 'a+') as save_file:
                for cec in cecs:
                    vec = list(np.array(cec) - np.array(com))
                    vec = [str(i) for i in vec]
                    save_file.write(' '.join(vec) + '\n')
                    save_file.flush()
        self.close_file()

    def compress_trj(self, line_range_list):
        line_range = list()
        for (i, j) in zip(line_range_list[0::2], line_range_list[1::2]):
            line_range += list(range(i-1, j))
        self.find_file(evb_flag=False)
        self.open_file(evb_flag=False)
        with open('%s%s.compressed.lammpstrj' % (self.work_dir, self.sys_name), 'w') as out_file:
            is_first_frame = True
            system_size = 0
            while True:
                self.trj_read_frame()
                if self.trj_time == 'EOF':
                    break
                if is_first_frame:
                    is_first_frame = False
                    system_size = len(self.trj_info)
                if not is_first_frame and system_size != len(self.trj_info):
                    break
                if self.trj_time % 5000 == 0:
                    print('[TRJ] %s%s DONE FRAME %s' % (self.work_dir, self.sys_name, self.trj_time))
                out_file.write('ITEM: TIMESTEP\n%s\n' % self.trj_time)
                out_file.write('ITEM: NUMBER OF ATOMS\n%s\n' % len(line_range))
                out_file.write('ITEM: BOX BOUNDS pp pp pp\n%s %s\n%s %s\n%s %s\nITEM: ATOMS id type mol x y z\n'
                               % (self.bounds[0][0], self.bounds[0][1],
                                  self.bounds[1][0], self.bounds[1][1],
                                  self.bounds[2][0], self.bounds[2][1]))
                out_file.flush()
                atom_id = 1
                for i, line in enumerate(self.trj_info):
                    if i in line_range:
                        line = line.strip('\n').split()
                        out_file.write('%s %s\n' % (atom_id, ' '.join(line[1:])))
                        atom_id += 1
        self.close_file(evb_flag=False)

    def compress_evb(self, frequency):
        self.find_file(trj_flag=False)
        self.open_file(trj_flag=False)
        with open('%s%s.compressed.evb' % (self.work_dir, self.sys_name), 'w') as out_file:
            while True:
                self.evb_read_frame()
                if not self.evb_info:
                    break
                if self.evb_time % 5000 == 0:
                    print('[EVB] %s%s DONE FRAME %s' % (self.work_dir, self.sys_name, self.evb_time))
                if self.evb_time % frequency == 0:
                    for line in self.evb_cache:
                        out_file.write(line)
                    out_file.flush()
        self.close_file(trj_flag=False)
