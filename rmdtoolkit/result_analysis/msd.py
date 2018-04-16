# Python 3.6.1

import os

import numpy as np

from rmdtoolkit.basic.result import Result
from rmdtoolkit.basic.tool import log_error, law_of_cosines, string_is_true


class MSD(Result):
    def __init__(self):
        super().__init__()
        self.sample_range = 1000
        self.sample_freq = 5
        self.stop_frame = float('inf')

        self.normal_flag = True
        self.com_correction = False
        self.lateral_flag = False
        self.radial_flag = False

        self.target_atoms = list()
        self.msd_result_normal = list()
        self.msd_result_radial = list()
        self.msd_result_lateral = list()

        self.sys_names = list()  # override auto searched results
        self.atom_tot = 0
        self.frame_tot = 0

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'SYSNAMES':  # should be in order
                            self.sys_names = line[1:]
                        elif line[0] == 'SAMPLE_RANGE':  # maximum x axis coordinate (delta t) in MSD figure
                            self.sample_range = int(line[1])
                        elif line[0] == 'SAMPLE_FREQ':  # sampling freq, delta delta_t while looping over delta_ts
                            self.sample_freq = int(line[1])
                        elif line[0] == 'NORMAL_FLAG':
                            self.normal_flag = True if string_is_true(line[1]) else False
                        elif line[0] == 'COM_FLAG':
                            self.com_correction = True if string_is_true(line[1]) else False
                        elif line[0] == 'RADIAL_FLAG':
                            self.radial_flag = True if string_is_true(line[1]) else False
                        elif line[0] == 'LATERAL_FLAG':
                            self.lateral_flag = True if string_is_true(line[1]) else False
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def calculate_msd(self):
        if not self.sys_names:
            self.sys_names = sorted([file_name.strip('.lammpstrj')
                                     for file_name in os.listdir(self.work_dir) if file_name.endswith('.lammpstrj')])
        atom_coords_list = list()
        coms_list = list()
        is_first_frame = True
        system_size = 0
        last_frame = 0
        for sys_name in self.sys_names:
            self.sys_name = sys_name
            self.find_file(evb_flag=self.evb_required)
            self.open_file(evb_flag=self.evb_required)
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
                if self.evb_required:
                    self.evb_match_trj()
                    if not self.evb_info:
                        break
                # if self.frame_tot > 5:  # For debug
                #     break
                last_frame = self.trj_time
                self.frame_tot += 1

                if self.frame_tot % 10 == 0:
                    print('[MSD] {} READING FRAME {}'.format(self.save_tag, self.trj_time))

                if is_first_frame:
                    is_first_frame = False
                    system_size = len(self.trj_info)

                temp_coords = list()
                for target_atom in self.target_atoms:
                    atom_coords = self.atom_find(target_atom)
                    temp_coords.extend(atom_coords)
                self.atom_tot = len(temp_coords)
                atom_coords_list.append(temp_coords)
                if self.com_correction or self.lateral_flag or self.radial_flag:
                    coms_list.append(self.find_com())
            self.close_file(evb_flag=self.evb_required)

        for delta_t in range(1, self.sample_range, self.sample_freq):
            if delta_t % 100 == 0:
                print('{} PROCESSING DELTA_T={}'.format(self.save_tag, delta_t))
            # if delta_t >= self.frame_tot - 1:
            #     break
            average_norm = 0.0
            average_rad = 0.0
            average_lat = 0.0
            # Loop over starting times.
            for t in range(self.frame_tot - delta_t):
                for i_atom in range(self.atom_tot):
                    if self.normal_flag:
                        if self.com_correction:
                            delta_r = (atom_coords_list[t][i_atom] - coms_list[t]) \
                                    - (atom_coords_list[t + delta_t][i_atom] - coms_list[t + delta_t])
                        else:
                            delta_r = atom_coords_list[t][i_atom] - atom_coords_list[t + delta_t][i_atom]
                        delta_r = np.squeeze(np.asarray(delta_r))
                        average_norm += np.dot(delta_r, delta_r)
                    if self.radial_flag:
                        delta_r1 = atom_coords_list[t][i_atom] - coms_list[t]
                        delta_r1 = np.squeeze(np.asarray(delta_r1))
                        delta_r1 = np.sqrt(np.dot(delta_r1, delta_r1))
                        delta_r2 = atom_coords_list[t + delta_t][i_atom] - coms_list[t + delta_t]
                        delta_r2 = np.squeeze(np.asarray(delta_r2))
                        delta_r2 = np.sqrt(np.dot(delta_r2, delta_r2))
                        average_rad += (delta_r1 - delta_r2) ** 2
                    if self.lateral_flag:
                        oa = atom_coords_list[t][i_atom] - coms_list[t]
                        ob = atom_coords_list[t + delta_t][i_atom] - coms_list[t + delta_t]
                        theta = law_of_cosines(oa, ob)
                        average_lat += theta
            # Normalize.
            try:
                if self.normal_flag:
                    average_norm /= self.atom_tot * (self.frame_tot - delta_t)
                    self.msd_result_normal.append(average_norm)
                if self.radial_flag:
                    average_rad /= self.atom_tot * (self.frame_tot - delta_t)
                    self.msd_result_radial.append(average_rad)
                if self.lateral_flag:
                    average_lat /= self.atom_tot * (self.frame_tot - delta_t)
                    self.msd_result_lateral.append(average_lat)
            except ZeroDivisionError:
                log_error(__file__, '{} {} {} {}'.format(average_norm, self.atom_tot, self.frame_tot, delta_t))

    def save_msd_result(self):
        d_list = range(1, self.sample_range, self.sample_freq)
        if self.normal_flag:
            with open('{}{}.msd'.format(self.save_dir, self.save_tag), 'w') as save_file:
                for d, r in zip(d_list, self.msd_result_normal):
                    save_file.write('{} {}\n'.format(d, r))
        if self.radial_flag:
            with open('{}{}.msd-rad'.format(self.save_dir, self.save_tag), 'w') as save_file:
                for d, r in zip(d_list, self.msd_result_radial):
                    save_file.write('{} {}\n'.format(d, r))
        if self.lateral_flag:
            with open('{}{}.msd-lat'.format(self.save_dir, self.save_tag), 'w') as save_file:
                for d, r in zip(d_list, self.msd_result_lateral):
                    save_file.write('{} {}\n'.format(d, r))
