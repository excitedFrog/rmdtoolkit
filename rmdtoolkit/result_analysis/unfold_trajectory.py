# Python 3.6.1

import os

import numpy as np

from basic.result import Result


class UnfoldTrj(Result):
    def __init__(self):
        super().__init__()
        self.sys_names = list()
        self.save_tag = str()

    def unfold_trajectory(self):
        os.chdir(self.work_dir)
        temp_coords = np.array([])
        new_coords = np.array([])
        f_matrix = np.array([])
        is_first_frame = True
        last_frame = 0
        with open('%s%sunfolded.lammpstrj' % (self.work_dir, self.save_tag), 'w') as out_file:
            for sys_name in self.sys_names:
                self.sys_name = sys_name
                self.find_file(evb_flag=False)
                self.open_file(evb_flag=False)
                while True:
                    self.trj_read_frame()
                    # Break and go to next file upon 'EOF' event.
                    if self.trj_time == 'EOF':
                        break
                    # Skip overlapping frames.
                    elif self.trj_time < last_frame:
                        continue
                    last_frame = self.trj_time

                    coords = np.array(
                        [np.array([float(line.split()[i]) for i in range(3, 6)]) for line in self.trj_info])
                    if is_first_frame:
                        is_first_frame = False
                        f_matrix = np.zeros(coords.shape)  # "f" stands for "Fold", indicating index of box
                        new_coords = coords
                    else:
                        d_matrix = coords - temp_coords  # "d" stands for "Displacement"
                        for i, row in enumerate(d_matrix):
                            for j, d in enumerate(row):
                                if abs(d) > 0.5 * self.box_len[j]:
                                    f_matrix[i][j] -= np.sign(d)
                                new_coords[i][j] = coords[i][j] + f_matrix[i][j] * self.box_len[j]
                    temp_coords = coords

                    # Write
                    out_file.write('ITEM: TIMESTEP\n%s\nITEM: NUMBER OF ATOMS\n%s\nITEM: BOX BOUNDS pp pp pp\n'
                                   % (self.trj_time, self.total_atoms))
                    for i in range(3):
                        out_file.write('%s %s\n' % (self.bounds[i][0], self.bounds[i][1]))
                    out_file.write('ITEM: ATOMS id type mol x y z\n')
                    for i, row in enumerate(new_coords):
                        line = self.trj_info[i].split()
                        out_file.write('%s %s %s %s %s %s\n' % (line[0], line[1], line[2], row[0], row[1], row[2]))
                    if self.trj_time % 10000 == 0:
                        print('Done Frame %s' % self.trj_time)
