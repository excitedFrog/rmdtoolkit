# Python 3.6.1

import numpy as np
import pandas as pd

from rmdtoolkit.basic.result import Result
from rmdtoolkit.basic.tool import pbc_dist


class Density(Result):
    def __init__(self):
        super().__init__()

        # settings
        self.directions = list()
        self.bin_size = 1.0
        self.sample_frequency = 1.0
        self.start_point = np.array([])
        self.stop_point = np.array([])
        self.target_atoms = list()
        self.stop_at = float('inf')
        # results
        self.box_volume = 0.
        self.frame_tot = 0
        self.sample_points = np.array([])
        self.density = np.array([])

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'DIRECTIONS':  # The degree of freedoms
                            self.directions.extend(line[1:])
                        elif line[0] == 'BIN_SIZE':  # Size of the bin for histogramic density
                            self.bin_size = float(line[1])
                        elif line[0] == 'SAMPLE_FREQ':
                            self.sample_frequency = float(line[1])
                        elif line[0] == 'START_POINT':
                            self.start_point = np.array(list(map(float, line[1:])))
                        elif line[0] == 'STOP_POINT':
                            self.stop_point = np.array(list(map(float, line[1:])))
                        elif line[0] == 'STOP_AT':
                            self.stop_at = int(line[1])
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))
        self.sample_points = np.vstack(map(np.ravel, np.meshgrid(*[np.arange(start, stop, self.sample_frequency)
                                       for start, stop in zip(self.start_point, self.stop_point)]))).T
        self.density = np.zeros(len(self.sample_points))

    def calculate_density(self, com_correction=False, method='histogram'):
        is_first_frame = True
        system_size = 0
        last_frame = 0
        self.find_file()
        self.open_file()
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            elif self.trj_time < last_frame:
                continue
            if len(self.trj_info) != system_size and not is_first_frame:
                continue
            if self.frame_tot > self.stop_at:
                break
            last_frame = self.trj_time
            if self.evb_flag:
                flag = self.evb_match_trj()
                if flag:
                    continue
            self.frame_tot += 1
            if self.frame_tot % 1 == 0:
                print('[DENSITY] Processing Frame {}'.format(self.trj_time))
            self.parse_trj_info()
            self.histogram_density(com_correction=com_correction)
        self.density /= self.frame_tot
        self.save_density()

    def calculate_density_radial(self, com_correction=False, method='histogram'):
        is_first_frame = True
        system_size = 0
        last_frame = 0
        self.find_file()
        self.open_file()
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            elif self.trj_time < last_frame:
                continue
            if len(self.trj_info) != system_size and not is_first_frame:
                continue
            if self.frame_tot > self.stop_at:
                break
            last_frame = self.trj_time
            if self.evb_flag:
                flag = self.evb_match_trj()
                if flag:
                    continue
            self.frame_tot += 1
            if self.frame_tot % 1 == 0:
                print('[DENSITY] Processing Frame {}'.format(self.trj_time))
            self.parse_trj_info()
            self.histogram_density_radial()

    def calculate_box_volume(self):
        xyz = [self.bin_size if __ in self.directions else self.box_len[_] for _, __ in enumerate(['x', 'y', 'z'])]
        self.box_volume = np.prod(xyz)

    def histogram_density(self, com_correction=False):
        self.calculate_box_volume()
        positions = self.atoms_find(df=True)
        positions = pd.DataFrame.as_matrix(positions, columns=self.directions)
        print(positions)
        exit()
        if com_correction:
            positions -= self.find_com()

        box_len = np.array([self.box_len[_] for _, __ in enumerate(['x', 'y', 'z']) if __ in self.directions])

        def number_in_box(point):
            select = [(pbc_dist(point, _, box_len, retarray=True) < self.bin_size).all() for _ in positions]
            return np.sum(select)

        temp_density = np.array(list(map(number_in_box, self.sample_points)), dtype=float)
        temp_density /= self.box_volume
        self.density += temp_density

    def histogram_density_radial(self):
        positions = self.atoms_find(df=True)
        positions = pd.DataFrame.as_matrix(positions, columns=self.directions)
        positions -= self.find_com()
        dists = np.linalg.norm(positions, axis=1).astype(str)

        with open('{}{}.density'.format(self.save_dir, self.save_tag), 'a') as save_file:
            save_file.write(' '.join(dists) + ' ')

    def kernel_density_gaussian(self):
        pass  # TODO

    def nearest_neighbor_density(self):
        pass  # TODO

    def save_density(self):
        with open('{}{}.density'.format(self.save_dir, self.save_tag), 'w') as save_file:
            save_file.write('{} density\n'.format(' '.join(self.directions)))
            for sample_point, density in zip(self.sample_points, self.density):
                save_file.write('{} {}\n'.format(' '.join(list(map(str, sample_point))), density))
