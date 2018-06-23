# Python 3.6.1

import os
import numpy as np
import pandas as pd
import itertools
from io import StringIO

from collections import OrderedDict

from rmdtoolkit.database.my_database import MyDB
from rmdtoolkit.basic.tool import string_is_true

# Pandas SettingWithCopyWarning is annoying, when you know what you are doing.
# I am turning it off here.
pd.options.mode.chained_assignment = None  # default='warn'


class Result(object):
    def __init__(self):
        super().__init__()

        self.trj_flag = True
        self.evb_flag = False
        self.input_flag = True

        self._sys_name = str()
        self.trj_basename = str()
        self.evb_basename = str()
        self.trj_extname = 'lammpstrj'
        self.evb_extname = 'evb'
        self.input_path = str()
        self.work_dir = str()

        self.tell_process_freq = 20
        self.stop_frame = float('inf')
        self.is_first_frame = True
        self.system_size = int()
        self.last_frame = int()
        self.frame_tot = int()

        self.trj_path = None
        self.trj_file = None
        self.trj_time = None
        self.trj_info = None
        self.parser_yes = False
        self.trj_info_parsed = dict()
        self.trj_info_df = pd.DataFrame()
        self.total_atoms = 0
        self.kw_list = list()  # kw means keyword
        self.kw_type = {'id': int, 'type': int, 'mol': int,
                        'x': float, 'y': float, 'z': float,
                        'vx': float, 'vy': float, 'vz': float,
                        'q': float}
        self.bounds = np.zeros((3, 2))
        self.box_len = np.zeros(3)
        self.system_volume = float()

        self.evb_path = None
        self.evb_file = None
        self.evb_time = None
        self.evb_info = None
        self.evb_cache = None

        self.atom_mass = dict()
        self.atom_tag = dict()
        self.total_species = 0
        self.specie_tags = list()
        self.specie_numbers = list()
        self.sections = OrderedDict()
        self.target_atoms = list()

        self.wat_o_type = int()
        self.wat_h_type = int()
        self.hyd_o_type = int()
        self.hyd_h_type = int()

        self.acceptor_center_type = int()
        self.donor_center_type = int()

        self._com_range_list = list()
        self.com_range = list()

        self.save_tag = 'default'
        self.save_dir = './'

        self.db = MyDB()

    @property
    def com_range_list(self):
        return self._com_range_list

    @com_range_list.setter
    def com_range_list(self, a_list):
        self._com_range_list = [int(_) for _ in a_list]
        for (i, j) in zip(self.com_range_list[0::2], self.com_range_list[1::2]):
            self.com_range += list(range(i, j+1))

    @property
    def sys_name(self):
        return self._sys_name

    @sys_name.setter
    def sys_name(self, string):
        self._sys_name = string
        self.trj_basename = self._sys_name
        self.evb_basename = self._sys_name

    def read_input(self):
        if not self.input_flag:
            return 0
        if not self.input_path:
            raise Exception('[ERROR] Input file not assigned!')
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'RMD':
                            self.evb_flag = True if string_is_true(line[1]) else False
                        elif line[0] == 'COMRANGE':
                            self.com_range_list = [int(_) for _ in line[1:]]
                        elif line[0] == 'SAVETAG':
                            self.save_tag = line[1]
                        elif line[0] == 'SAVEDIR':
                            self.save_dir = line[1]
                        elif line[0] == 'TPFREQ':
                            self.tell_process_freq = int(line[1]) if int(line[1]) > 0 else 20
                        elif line[0] == 'STOP_FRAME':  # stop frame (unit is frame, not simulation time)
                            self.stop_frame = int(line[1])

                        elif line[0] == 'SYS_NAME':
                            self.sys_name = line[1]
                        elif line[0] == 'TRJ_NAME':
                            self.trj_basename = line[1]
                        elif line[0] == 'EVB_NAME':
                            self.evb_basename = line[1]
                        elif line[0] == 'TRJ_EXTNAME':
                            self.trj_extname = line[1]
                        elif line[0] == 'EVB_EXTNAME':
                            self.evb_extname = line[1]

                        elif line[0] == 'SPECIE':
                            self.specie_tags.append(line[1])
                            self.specie_numbers.append(int(line[2]))
                            self.total_species += 1
                        elif line[0] == 'ATOM':  # Target atoms to be found
                            self.target_atoms.append(line[1:])
                        elif line[0] == 'WAT_O_TYPE':
                            self.wat_o_type = int(line[1])
                        elif line[0] == 'WAT_H_TYPE':
                            self.wat_h_type = int(line[1])
                        elif line[0] == 'HYD_O_TYPE':
                            self.hyd_o_type = int(line[1])
                        elif line[0] == 'HYD_H_TYPE':
                            self.hyd_h_type = int(line[1])
                        elif line[0] == 'ACP_CENT_TYPE':
                            self.acceptor_center_type = int(line[1])
                        elif line[0] == 'DON_CENT_TYPE':
                            self.donor_center_type = int(line[1])
                        elif line[0] == 'ATOM_SYMBOL':  # ATOM_SYMBOL type1, symbol1, type2, symbol2, ...
                            for atom_type, atom_tag in zip(line[1::2], line[2::2]):
                                self.atom_tag[int(atom_type)] = atom_tag
                                self.db.cursor.execute('SELECT ElementMass FROM Elements WHERE ElementName=\'{}\''
                                                       .format(atom_tag))
                                atom_mass = self.db.cursor.fetchone()[0]
                                self.atom_mass[int(atom_type)] = atom_mass
                        # atom mass is generated upon reading in atom symbol
                        # Write ATOM_MASS line after ATOM_SYMBOL line to override masses of selected atoms
                        elif line[0] == 'ATOM_MASS':  # ATOM_MASS type1, mass1, type2, mass2, ...
                            for atom_type, atom_mass in zip(line[1::2], line[2::2]):
                                self.atom_mass[int(atom_type)] = float(atom_mass)
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def find_file(self):
        if self.trj_flag:
            self.trj_path = '{}{}.{}'.format(self.work_dir, self.trj_basename, self.trj_extname)
            if not os.path.isfile(self.trj_path):
                raise Exception('File {} does not exist!'.format(self.trj_path))
        if self.evb_flag:
            self.evb_path = '{}{}.{}'.format(self.work_dir, self.evb_basename, self.evb_extname)
            if not os.path.isfile(self.evb_path):
                raise Exception('File {} does not exist!'.format(self.evb_path))

    def open_file(self):
        if self.trj_flag:
            self.trj_file = open(self.trj_path, 'r')
        if self.evb_flag:
            self.evb_file = open(self.evb_path, 'r')

    def close_file(self):
        if self.trj_flag:
            self.trj_file.close()
        if self.evb_flag:
            self.evb_file.close()

    def checked_read(self, parse=True):
        if not any((self.evb_flag, self.trj_flag)):
            raise Exception('At least one of evb_flag or trj_flag should be set to True.')
        if self.trj_flag:
            self.trj_read_frame()
            if not self.is_first_frame and len(self.trj_info) != self.system_size:
                return -1
            if self.trj_time == 'EOF':
                print('EOF')
                return -1
            if self.frame_tot > self.stop_frame:
                return -1
            elif self.trj_time < self.last_frame:
                return 1
            if self.evb_flag:
                self.evb_match_trj()
                if not self.evb_info:
                    return -1
            self.last_frame = self.trj_time
            if self.is_first_frame:
                self.is_first_frame = False
                self.system_size = len(self.trj_info)
            if parse:
                self.parse_trj_info()
            self.frame_tot += 1
            return 0
        else:
            self.evb_read_frame()
            if self.frame_tot > self.stop_frame:
                return -1
            elif self.evb_time < self.last_frame:
                return 1
            if not self.evb_info:
                return -1
            self.frame_tot += 1
            return 0

    def tell_process(self):
        if self.frame_tot % self.tell_process_freq == 0:
            print('[{}] {} READING FRAME {}'.format(self.__class__.__name__, self.save_tag, self.frame_tot))

    def void_func(self):
        pass

    def analysis_template(self, inner_compute_func, inner_save_func, outer_compute_func, outer_save_func,
                          tell_process=True, parse=True):
        self.read_input()
        self.find_file()
        self.open_file()
        while True:
            checksum = self.checked_read(parse=parse)
            if checksum == 1:
                continue
            elif checksum == -1:
                break
            if tell_process:
                self.tell_process()
            inner_compute_func()
            inner_save_func()
        outer_compute_func()
        outer_save_func()
        self.close_file()

    # =================================================================================================================
    # || InfoGet Methods ||
    # =================================================================================================================
    def trj_read_frame(self):
        trj_info = list()
        trj_time = 'EOF'
        break_trigger = 0
        while True:
            last_pos = self.trj_file.tell()
            line = self.trj_file.readline()
            if not line:
                break
            lines = line.split()
            if lines[0] == 'ITEM:':
                if lines[1] == 'TIMESTEP':
                    break_trigger += 1
                    if break_trigger == 1:
                        trj_time = int(self.trj_file.readline().strip('\n'))
                    elif break_trigger > 1:
                        self.trj_file.seek(last_pos)
                        break
                elif lines[1] == 'BOX':
                    for i in range(3):
                        bounds = self.trj_file.readline().split()
                        self.bounds[i][0] = float(bounds[0])
                        self.bounds[i][1] = float(bounds[1])
                        self.box_len[i] = float(bounds[1]) - float(bounds[0])
                    self.system_volume = self.box_len[0] * self.box_len[1] * self.box_len[2]
                elif lines[1] == 'NUMBER':
                    self.total_atoms = int(self.trj_file.readline().strip('\n'))
                elif lines[1] == 'ATOMS':
                    if not self.parser_yes:
                        self.parser_yes = True
                        for item in lines[2:]:
                            self.trj_info_parsed[item] = list()
                            self.kw_list.append(item)
            else:
                trj_info.append(line.strip('\n'))
        self.trj_info = trj_info
        self.trj_time = trj_time
        self.parse_trj_info()

    def evb_read_frame(self):
        self.evb_cache = list()
        evb_info = {'Time': 0,
                    'ComplexCount': 0,
                    'ReactionCenters': list(),
                    'SHELLS': list(),
                    'EigenVectors': list(),
                    'CECS': list()}
        evb_time = 'EOF'
        complex_count = ''
        while True:
            line = self.evb_file.readline()
            if not line:
                break
            if line.strip('\n') == 'END_OF_COMPLEX %s' % str(complex_count):
                self.evb_cache.append(line)
                break
            else:
                self.evb_cache.append(line)
                line = line.split()
                if len(line) > 0:
                    # Get Timestep
                    if line[0] == 'TIMESTEP':
                        try:
                            evb_time = int(line[1].strip('\n'))
                        except IndexError:
                            break
                        evb_info['Time'] = evb_time
                    # Get Complex Count
                    elif line[0] == 'COMPLEX_COUNT':
                        try:
                            complex_count = int(line[1])
                        except IndexError:
                            break
                        evb_info['ComplexCount'] = complex_count
                    # Get Reaction Centers
                    elif line[0] == 'REACTION_CENTER_LOCATION':
                        while True:
                            line = self.evb_file.readline()
                            if not line:
                                break
                            self.evb_cache.append(line)
                            line = line.split()
                            evb_info['ReactionCenters'].append([int(item) for item in line])
                            if int(line[0]) == complex_count:
                                break
                    # Get Eigen Vectors
                    elif line[0] == 'EIGEN_VECTOR':
                        line = self.evb_file.readline()
                        self.evb_cache.append(line)
                        line = line.split()
                        if not line:
                            break
                        eigen_vector = [float(item) for item in line]
                        eigen_vector.sort(reverse=True)
                        evb_info['EigenVectors'].append(eigen_vector)
                    # Get CEC Coordinates
                    elif line[0] == 'CEC_COORDINATE':
                        line = self.evb_file.readline()
                        if not line:
                            break
                        self.evb_cache.append(line)
                        line = line.split()
                        evb_info['CECS'].append([float(item) for item in line])
                    # Get Solvation Shell Info
                    elif line[0] == 'STATES':
                        shells = list()
                        while True:
                            line = self.evb_file.readline()
                            self.evb_cache.append(line)
                            line = line.strip('\n').split()
                            if not line:
                                break
                            if not line[0].isdigit():
                                break
                            shells.append([int(item) for item in line])
                        evb_info['SHELLS'].append(shells)
        # Check evb information integrity.
        info_list = [evb_info[_] for _ in evb_info]
        if all(info_list) or evb_time == 0:
            evb_info['ReactionCenters'] = np.array(evb_info['ReactionCenters'])
            evb_info['EigenVectors'] = np.array(evb_info['EigenVectors'])
            evb_info['CECS'] = np.array(evb_info['CECS'])
            evb_info['SHELLS'] = np.array(evb_info['SHELLS'])
            self.evb_info = evb_info
            self.evb_time = evb_time
        else:
            self.evb_info = None  # break by call "if not self.evb_info: break"

    def evb_match_trj(self):
        while True:
            self.evb_read_frame()
            if self.evb_time == self.trj_time:
                break
            if not self.evb_time:
                continue
            if self.evb_time > self.trj_time:
                return 1

    def parse_trj_info(self):
        info = [' '.join(self.kw_list)]
        info.extend(self.trj_info)
        self.trj_info_df = pd.read_csv(StringIO('\n'.join(info)), delim_whitespace=True, dtype=self.kw_type)

    @staticmethod
    def get_pos(df):
        return pd.DataFrame.as_matrix(df, columns=['x', 'y', 'z'])

    @staticmethod
    def get_pos_df(df):
        return pd.DataFrame(df, columns=['x', 'y', 'z'])

    @staticmethod
    def get_vel(df):
        return pd.DataFrame.as_matrix(df, columns=['vx', 'vy', 'vz'])

    @staticmethod
    def get_vel_df(df):
        return pd.DataFrame(df, columns=['vx', 'vy', 'vz'])

    def get_mass(self, df):
        types = df.as_matrix(columns=['type'])
        return np.array([self.atom_mass[_[0]] for _ in types])

    # =================================================================================================================
    # ||AtomFind Methods||
    # =================================================================================================================
    def atoms_find(self, target_atoms=None, df=False):  # Target atoms is a list of list.
        temp = list()
        if not target_atoms:
            target_atoms = self.target_atoms

        for target_atom in target_atoms:
            temp.append(self.atom_find(target_atom, df=df))
        if df:
            return pd.concat(temp)
        else:
            return np.concatenate(temp)

    def atom_find(self, target_atom, df=False):
        func_dict = {'COM': self.find_com,
                     'DON': self.find_donor,
                     'CEC': self.find_cec}
        if target_atom[0] in func_dict:
            return func_dict[target_atom[0]](df=df)
        else:
            mol_tag = target_atom[0]
            atom_ids = [int(i) for i in target_atom[1:-1]]
            return self.find_atom(mol_tag, atom_ids, df=df)

    def find_com(self, df=False):
        info = self.trj_info_df[self.trj_info_df['id'].isin(self.com_range)] if self.com_range else self.trj_info_df
        masses = self.get_mass(info)
        positions = self.get_pos(info)
        com = np.dot(masses, positions) / np.sum(masses)
        if df:
            return pd.DataFrame(np.array([com]), columns=['x', 'y', 'z'])
        else:
            return np.array(com)

    def find_donor(self, df=False):
        info = self.trj_info_df[self.trj_info_df['mol'].isin(self.evb_info['ReactionCenters'][:, 1])
                                & self.trj_info_df['type'].isin([self.wat_o_type, self.hyd_o_type])]
        if df:
            return self.get_pos_df(info)
        else:
            return self.get_pos(info)

    def find_cec(self, df=False):
        if df:
            return pd.DataFrame(np.array(self.evb_info['CECS']), columns=['x', 'y', 'z'])
        else:
            return np.array(self.evb_info['CECS'])

    # find_atom takes type now, instead of atom_id
    def find_atom(self, mol_tag, atom_types, df=False):
        if not self.sections:
            self.generate_sections()
        if mol_tag == 'WAT' and 'HYD' in self.specie_tags:
            info = self.trj_info_df[self.trj_info_df['mol'].isin(itertools.chain(self.sections['WAT'],
                                                                                 self.sections['HYD']))
                                    & self.trj_info_df['type'].isin(atom_types)]
        else:
            info = self.trj_info_df[self.trj_info_df['mol'].isin(self.sections[mol_tag])
                                    & self.trj_info_df['type'].isin(atom_types)]
        if df:
            return self.get_pos_df(info)
        else:
            return self.get_pos(info)

    # =================================================================================================================
    # ||Misc||
    # =================================================================================================================
    def fetch_mol(self, mol_ids):
        mol_lines = list()
        for line in self.trj_info:
            if int(line.split()[2]) in mol_ids:
                mol_lines.append(line)
        return mol_lines

    # Note: This now uses the mol_id, instead of line_number, so the full specie setting should be included in input,
    # instead of just printed species.
    def generate_sections(self):
        for i, (specie, number) in enumerate(zip(self.specie_tags, self.specie_numbers)):
            if i == 0:
                self.sections[specie] = list(range(1, number+1))
            else:
                start = self.sections[self.specie_tags[i-1]][-1]
                self.sections[specie] = list(range(start+1, start+number+1))

    def lammps_trj_frame_to_pdb(self):
        out_path = '%s%snew.pdb' % (self.work_dir, self.sys_name)
        ref_path = '%s%s.pdb' % (self.work_dir, self.sys_name)
        out_lines = ['HEADER\nTITLE NONE\nREMARK NONE\n']

        ref_lines = list()
        with open(ref_path, 'r') as ref_file:
            for line in ref_file:
                if line.startswith('ATOM'):
                    ref_lines.append(line)

        for trj_line, ref_line in zip(self.trj_info, ref_lines):
            trj_l = trj_line.split()
            ref_l = ref_line.split()
            serial = ref_l[1]
            name = ref_l[2]
            res_name = ref_l[3]
            chain_id = ref_l[4]
            res_seq = ref_l[5]
            x, y, z = [round(float(num), 3) for num in trj_l[3:6]]
            occupancy = ref_l[9]
            temp_factor = ref_l[10]
            if len(ref_l[11]) > 2:
                element = ref_l[11][:-2]
                charge = ref_l[11][-2:]
            else:
                element = ref_l[11]
                charge = ''

            pdb_line = 'ATOM{serial}{name}{altLoc}{resName}{chainID}{resSeq}{iCode}' \
                       '{x}{y}{z}{occupancy}{tempFactor}{element}{charge}\n'\
                .format(serial=serial.rjust(7), name=name.rjust(5), altLoc=''.rjust(1), resName=res_name.rjust(3),
                        chainID=chain_id.rjust(2), resSeq=res_seq.rjust(4), iCode=''.rjust(1), x=str(x).rjust(11),
                        y=str(y).rjust(8), z=str(z).rjust(8), occupancy=occupancy.rjust(6),
                        tempFactor=temp_factor.rjust(6), element=element.rjust(12), charge=charge.rjust(2))
            out_lines.append(pdb_line)

        out_lines.append('END')
        with open(out_path, 'w') as out_file:
            for out_line in out_lines:
                out_file.write(out_line)
