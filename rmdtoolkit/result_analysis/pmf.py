# Python 3.6.1

import os

import numpy as np

from basic.result import Result
from basic.tool import log_error
from result_analysis.water_cluster import WaterCluster


class PMF(WaterCluster, Result):
    def __init__(self):
        super().__init__()
        self.sample_freq = 1
        self.bin_size = 0.001
        self.delta_range = 0.6

        self.frame_tot = 0
        self.complex_tot = 0
        self.n_bin = 0
        self.counts = None

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'BIN_SIZE':
                            self.bin_size = float(line[1])
                        elif line[0] == 'XI_SAMPLE_FREQ':
                            self.sample_freq = int(line[1])
                        elif line[0] == 'DELTA_RANGE':
                            self.delta_range = float(line[1])
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def count_xi(self):
        self.find_file(trj_flag=False)
        self.open_file(trj_flag=False)

        is_first_frame = True
        self.n_bin = int(1.0/self.bin_size)
        while True:
            self.evb_read_frame()
            if not self.evb_info or self.evb_time == 'EOF':
                break
            if self.evb_time % 100 == 0:
                print('{} READING FRAME {}'.format(self.save_tag, self.evb_time))

            eigen_vectors = self.evb_info['EigenVectors']
            if is_first_frame:
                is_first_frame = False
                self.complex_tot = len(eigen_vectors)
                self.counts = np.zeros((self.complex_tot, self.n_bin))

            # Sample frequency
            if not self.evb_time % self.sample_freq == 0:
                continue

            self.frame_tot += 1
            for i_complex, eigen_vector in enumerate(eigen_vectors):
                try:
                    xi = eigen_vector[0] ** 2 - eigen_vector[1] ** 2
                except IndexError:  # Naked hydronium
                    xi = 1.0 - self.bin_size
                i_bin = int(xi / self.bin_size)
                try:
                    self.counts[i_complex][i_bin] += 1
                except KeyError:
                    log_error(__file__, '{} {}\n'.format(i_complex, i_bin))
                    self.frame_tot -= 1
                except IndexError:
                    log_error(__file__, '{} {}\n'.format(i_complex, i_bin))
                    self.frame_tot -= 1
        self.close_file(trj_flag=False)

    def count_delta(self):
        os.chdir(self.work_dir)
        self.find_file()
        self.open_file()
        self.n_bin = int(self.delta_range / self.bin_size)

        system_size = 0
        is_first_frame = True
        while True:
            self.trj_read_frame()
            if not is_first_frame and len(self.trj_info) != system_size:
                break
            if self.trj_time == 'EOF':
                break
            self.evb_match_trj()
            if not self.evb_info:
                break
            if is_first_frame:
                is_first_frame = False
                system_size = len(self.trj_info)
                self.complex_tot = self.evb_info['ComplexCount']
                self.counts = np.zeros((self.complex_tot, self.n_bin))

            self.frame_tot += 1
            if self.frame_tot % 10 == 0:
                print('{} READING FRAME {}'.format(self.save_tag, self.trj_time))

            clusters, center_ids = self.get_evb_clusters()
            for i_complex, (cluster, center_id) in enumerate(zip(clusters, center_ids)):
                delta = self.calculate_psp(cluster, center_id)
                i_bin = int(delta / self.bin_size)
                try:
                    self.counts[i_complex][i_bin] += 1
                except IndexError:
                    log_error(__file__, '{} {}\n'.format(i_complex, i_bin))
        self.close_file()

    def save_pmf_result(self):
        for i_complex in range(self.complex_tot):
            with open('%s%s-comp%s.pmf' % (self.save_dir, self.save_tag, str(i_complex+1)), 'w') as save_file:
                save_file.write('%s\n' % self.frame_tot)
                for i_bin in range(self.n_bin):
                    save_file.write('{:.5f} {}\n'.format(i_bin * self.bin_size, self.counts[i_complex][i_bin]))
