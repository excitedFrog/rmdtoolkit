# Python 3.6.1

import os

import numpy as np

from rmdtoolkit.basic.result import Result
from rmdtoolkit.basic.tool import pbc_dist


class RDF(Result):
    def __init__(self):
        super().__init__()
        self.rdf_input_file = None

        self.cent_atom = list()
        self.cent_name = None
        self.dstr_atoms = list()

        self.count = dict()
        self.bin_size = 0.01
        self.max_dist = 60.0
        self.n_bin = int(self.max_dist / self.bin_size)
        self.frame_tot = 0
        self.cent_tot = dict()
        self.normalise = dict()

        self.evb_required = False
        self.complex_flag = False
        self.complex_cutoff = 3.5

    def read_input(self):
        super().read_input()
        with open(self.rdf_input_file, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'CENT':
                            self.cent_atom = line[1:]
                            self.cent_name = self.cent_atom[-1]
                        elif line[0] == 'DSTR':
                            self.dstr_atoms.append(line[1:])
                            dstr_name = line[-1]
                            self.cent_tot[dstr_name] = 0
                            self.normalise[dstr_name] = 0
                        elif line[0] == 'BIN_SIZE':
                            self.bin_size = float(line[1])
                        elif line[0] == 'MAX_DIST':
                            self.max_dist = float(line[1])
                        # in inhomogeneous coordinated systems, turn on this feature, it allows you to set cutoff
                        # distance to judge if a center atom is part of a complex, which will change the normalization
                        # factor, giving back a correct coordination number
                        elif line[0] == 'COMPLEX_FLAG':
                            self.complex_flag = bool(line[1] in ['True', 'true', 'Yes', 'yes'])
                        elif line[0] == 'COMPLEX_CUTOFF':
                            self.complex_flag = True
                            self.complex_cutoff = float(line[1])
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))
        self.n_bin = int(self.max_dist / self.bin_size)
        for target_atom in self.cent_atom + self.dstr_atoms:
            if target_atom[0] in ['CEC', 'O*', 'H*', 'OFSS']:
                self.evb_required = True

    def count_rdf(self):
        self.find_file(evb_flag=self.evb_required)
        self.open_file(evb_flag=self.evb_required)

        for dstr_atom in self.dstr_atoms:
            dstr_name = dstr_atom[-1]
            self.count[dstr_name] = np.zeros(self.n_bin)

        is_first_frame = True
        system_size = 0
        last_frame = 0
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            elif self.trj_time < last_frame:
                continue
            if len(self.trj_info) != system_size and not is_first_frame:
                continue
            last_frame = self.trj_time
            if self.evb_required:
                flag = self.evb_match_trj()
                if flag:
                    continue
            self.frame_tot += 1

            if self.frame_tot % 10 == 0:
                print('[RDF] {} PROCESSING FRAME {}'.format(self.save_tag, self.trj_time))

            dstr_coords = dict()
            if is_first_frame:
                is_first_frame = False
                system_size = len(self.trj_info)

            center_coords = self.atom_find(self.cent_atom)
            for dstr_atom in self.dstr_atoms:
                dstr_name = dstr_atom[-1]
                dstr_coords[dstr_name] = self.atom_find(dstr_atom)

            for dstr_name in dstr_coords:
                self.cent_tot[dstr_name] = 0
                for cent_coord in center_coords:
                    is_complex = False
                    for dstr_coord in dstr_coords[dstr_name]:
                        r_real = pbc_dist(cent_coord, dstr_coord, self.box_len)
                        i_bin = int(r_real / self.bin_size)
                        if self.complex_flag:
                            if r_real < self.complex_cutoff:
                                is_complex = True
                                try:
                                    self.count[dstr_name][i_bin] += 1
                                except IndexError:
                                    pass
                        else:
                            try:
                                self.count[dstr_name][i_bin] += 1
                            except IndexError:
                                pass
                    if self.complex_flag:
                        if is_complex:
                            self.cent_tot[dstr_name] += 1
                    else:
                        self.cent_tot[dstr_name] += 1
                # print(self.save_tag, dstr_name, self.cent_tot[dstr_name])
                self.normalise[dstr_name] += self.cent_tot[dstr_name]
        self.close_file(evb_flag=self.evb_required)

    def save_rdf_result(self):
        for dstr_atom in self.dstr_atoms:
            dstr_name = dstr_atom[-1]
            with open('{}{}-{}-{}.rdf'.format(self.save_dir, self.save_tag, self.cent_name, dstr_name), 'w')\
                    as save_file:
                save_file.write('{}\n'.format(self.normalise[dstr_name]))
                for i, val in enumerate(self.count[dstr_name]):
                    save_file.write('{:.5f} {}\n'.format(i*self.bin_size, val))
