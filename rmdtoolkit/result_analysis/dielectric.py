# Python 3.6.1

import os
import itertools

import numpy as np

from basic.tool import chunks
from basic.result import Result


class Dielectric(Result):
    def __init__(self):
        super().__init__()
        self.sys_names = list()
        self.evb_required = False
        self.total_frames = 0

        self.times = list()
        self.component_types = list()  # ion or molecule
        self.component_mols = list()  # tag of molecule to be calculated. consistent with SPECIE tag
        self.dipole_results = dict()

        self.cross_pairs = list()  # cross correlation functions to be calculated
        self.cross_results = dict()
        self.max_dt = 2000000  # unit of dt is "frame"

    def work1(self):
        self.read_input()
        self.generate_sections()
        self.calculate_dipoles()
        self.save_dipoles()

    def work2(self):
        self.read_input()
        self.load_dipoles()
        self.calculate_cross()
        self.save_cross()

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'DIELECTRIC_COMPONENT':
                            if line[1] in ['MOL', 'ION']:
                                self.component_types.append(line[1])
                            else:
                                print('!UNRECOGNIZED MOLECULE TYPE!')
                                exit()
                            mol_tag = line[2]
                            self.component_mols.append(mol_tag)
                            self.dipole_results[mol_tag] = list()
                        elif line[0] == 'DIELECTRIC_CROSS':
                            self.cross_pairs.append((line[1], line[2]))
                            self.cross_results['{}-{}'.format(line[1], line[2])] = list()
                    except IndexError:
                        raise Exception('No/incomplete value(s) found after tag \'{}\' in input file.'.format(line[0]))

    def calculate_dipoles(self):
        os.chdir(self.work_dir)
        is_first_frame = True
        system_size = 0
        last_frame = 0
        self.find_file(evb_flag=self.evb_required)
        self.open_file(evb_flag=self.evb_required)
        while True:
            self.trj_read_frame()
            if not is_first_frame and len(self.trj_info) != system_size:
                break
            if self.trj_time == 'EOF':
                break
            elif self.trj_time < last_frame:
                continue
            if self.evb_required:
                self.evb_match_trj()
                if not self.evb_info:
                    break
            # if self.total_frames > 100:  # For debug
            #     break
            last_frame = self.trj_time
            self.total_frames += 1

            if self.total_frames % 100 == 0:
                print('{} AT CALCULATING DIPOLE OF FRAME {}'.format(self.save_tag, self.trj_time))

            self.times.append(self.trj_time)
            result = 'NULL'
            for comp_type, comp_mol in zip(self.component_types, self.component_mols):
                # construct charge list from database
                self.db.cursor.execute('SELECT MoleculeID FROM Molecules WHERE MoleculeType=\'{}\''
                                       .format(comp_mol))
                molecule_id = self.db.cursor.fetchone()[0]
                self.db.cursor.execute('SELECT AtomID FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
                molecule_len = len(self.db.cursor.fetchall())
                charge_list = list()
                for atom_id in range(1, molecule_len+1):
                    self.db.cursor.execute('SELECT Charge FROM Atoms WHERE MoleculeID={} AND AtomID={}'
                                           .format(molecule_id, atom_id))
                    charge_list.append(self.db.cursor.fetchone()[0])
                charge_list = np.array(charge_list)
                if comp_type == 'MOL':
                    result = self.dipole_mol(comp_mol, charge_list, molecule_len)
                elif comp_type == 'ION':
                    result = self.current_ion(comp_mol, charge_list, molecule_len)
                self.dipole_results[comp_mol].append(result)

    def save_dipoles(self):
        with open('{}{}.dipoles'.format(self.save_dir, self.save_tag), 'w') as save_file:
            for i, time in enumerate(self.times):
                save_line = [str(time)]
                for comp in self.component_mols:
                    save_line.append(str(self.dipole_results[comp][i]))
                save_line.append('\n')
                save_file.write(' '.join(save_line))

    def load_dipoles(self):
        file_names = os.listdir(self.save_dir)
        dipole_file_names = [_ for _ in file_names if self.save_tag in _ and _.endswith('.dipoles')]
        [self.load_dipoles_loader(self.save_dir + _) for _ in dipole_file_names]
        self.times, unique_indices = np.unique(self.times, return_index=True)
        self.dipole_results = {_: np.take(self.dipole_results[_], unique_indices, axis=0)
                               for _ in self.dipole_results.keys()}
        self.total_frames = len(self.times)

    def load_dipoles_loader(self, path):
        with open(path, 'r') as load_file:
            for line in load_file:
                line = line.replace('[', '').replace(']', '')
                line = line.split()
                line = list(map(float, line))
                self.times.append(int(line[0]))
                arrays = np.array(chunks(line[1:], 3))
                [self.dipole_results[_].append(__) for _, __ in zip(self.component_mols, arrays)]

    def calculate_cross(self):
        for pair in self.cross_pairs:
            for dt in range(0, self.max_dt):
                if dt % 1 == 0:
                    print('{} AT PROCESSING PAIR {} DT {}'.format(self.save_tag, pair, dt))
                l1 = np.array(self.dipole_results[pair[0]]) - np.average(self.dipole_results[pair[0]])
                l2 = np.array(self.dipole_results[pair[1]]) - np.average(self.dipole_results[pair[1]])
                value = np.average(l1[0:self.total_frames-dt] * l2[dt:self.total_frames])
                self.cross_results['{}-{}'.format(pair[0], pair[1])].append(value)

    def save_cross(self):
        with open('{}{}.cross'.format(self.save_dir, self.save_tag), 'w') as save_file:
            for dt in range(0, self.max_dt):
                save_line = [str(dt)]
                for pair in self.cross_pairs:
                    save_line.append(str(self.cross_results['{}-{}'.format(pair[0], pair[1])][dt]))
                save_line.append('\n')
                save_file.write(' '.join(save_line))

    def dipole_mol(self, comp, charge_list, molecule_len):
        if comp == 'WAT':  # because with RMD the sequence of water and hydronium atoms mutates
            hyd_ids = self.evb_info['ReactionCenters'][:, 1]
            info = self.trj_info_df[self.trj_info_df['mol'].isin(itertools.chain(self.sections['WAT'],
                                                                                 self.sections['HYD']))
                                    & ~self.trj_info_df['mol'].isin(hyd_ids)]
            o_pos = self.get_pos(info[info['type'].isin([self.wat_o_type, self.hyd_o_type])])
            h_pos = self.get_pos(info[info['type'].isin([self.wat_h_type, self.hyd_h_type])])
            dipole = (np.sum(o_pos * self.wat_o_charge, axis=0) + np.sum(h_pos * self.wat_h_charge, axis=0))
        else:
            info = self.trj_info_df[self.trj_info_df['mol'].isin(self.sections[comp])]
            calibrate = info['id'].iloc[0] - 1
            info['id'] = info['id'].apply(lambda x: x-calibrate)
            info['dipole_x'] = info['x'] * charge_list[info['id'] % molecule_len - 1]
            info['dipole_y'] = info['y'] * charge_list[info['id'] % molecule_len - 1]
            info['dipole_z'] = info['z'] * charge_list[info['id'] % molecule_len - 1]
            dipole = np.sum(info[['dipole_x', 'dipole_y', 'dipole_z']].as_matrix(), axis=0)
        return dipole

    def current_ion(self, comp, charge_list, molecule_len):
        if comp == 'HYD':  # because with RMD the sequence of water and hydronium atoms mutates
            hyd_ids = self.evb_info['ReactionCenters'][:, 1]
            info = self.trj_info_df[self.trj_info_df['mol'].isin(itertools.chain(self.sections['WAT'],
                                                                                 self.sections['HYD']))
                                    & self.trj_info_df['mol'].isin(hyd_ids)]
            o_vel = self.get_vel(info[info['type'].isin([self.wat_o_type, self.hyd_o_type])])
            h_vel = self.get_vel(info[info['type'].isin([self.wat_h_type, self.hyd_h_type])])
            current = (np.sum(o_vel * self.wat_o_charge, axis=0) + np.sum(h_vel * self.wat_h_charge, axis=0))
        else:
            info = self.trj_info_df[self.trj_info_df['mol'].isin(self.sections[comp])]
            calibrate = info['id'].iloc[0] - 1
            info['id'] = info['id'].apply(lambda x: x - calibrate)
            info['current_x'] = info['vx'] * charge_list[info['id'] % molecule_len - 1]
            info['current_y'] = info['vy'] * charge_list[info['id'] % molecule_len - 1]
            info['current_z'] = info['vz'] * charge_list[info['id'] % molecule_len - 1]
            current = np.sum(info[['current_x', 'current_y', 'current_z']].as_matrix(), axis=0)
        return current

