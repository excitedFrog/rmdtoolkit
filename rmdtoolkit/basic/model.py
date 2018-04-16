# Python 3.6.1

import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from collections import OrderedDict

from mendeleev import element as elem

from rmdtoolkit.database.my_database import MyDB


class Model(object):
    def __init__(self):
        super().__init__()
        self.input_path = str()
        self.total_species = 0
        self.specie_tags = list()
        self.specie_numbers = list()
        self.sections = OrderedDict()

        self.psf_path = None
        self._pdb_path = None
        self._mol2_path = None
        self._data_path = None
        self.sys_name = None
        self.work_dir = None

        self.leaprc_path = None

        self.pdb_kw_type = {'serial': int, 'name': str, 'altLoc': str, 'resName': str, 'chainID': str, 'resSeq': int,
                            'iCode': str, 'x': float, 'y': float, 'z': float, 'occupancy': str, 'tempFactor': str,
                            'element': str, 'charge': str}
        self.pdb_kw_range = {'serial': [6, 11], 'name': [13, 16], 'altLoc': [16, 17], 'resName': [17, 20],
                             'chainID': [21, 22], 'resSeq': [22, 26], 'iCode': [26, 27],
                             'x': [30, 38], 'y': [38, 46], 'z': [46, 54],
                             'occupancy': [54, 60], 'tempFactor': [60, 66], 'element': [76, 78], 'charge': [78, 80]}
        self.pdb_kw_list = ['serial', 'name', 'altLoc', 'resName', 'chainID', 'resSeq', 'iCode', 'x', 'y', 'z',
                            'occupancy', 'tempFactor', 'element', 'charge']

        self.pdb_df = pd.DataFrame()
        self.data_atom_df = pd.DataFrame()
        self.data_bond_df = pd.DataFrame()
        self.data_angle_df = pd.DataFrame()
        self.data_dihedral_df = pd.DataFrame()

        self._com_range_list = list()
        self.com_range = list()

        self.db = MyDB()

    @property
    def pdb_path(self):
        return self._pdb_path

    @pdb_path.setter
    def pdb_path(self, path):
        self._pdb_path = path
        self.sys_name = path.split('/')[-1].split('.')[0]
        self.work_dir = path[0:-len(path.split('/')[-1])]

    @property
    def mol2_path(self):
        return self._mol2_path

    @mol2_path.setter
    def mol2_path(self, path):
        self._mol2_path = path
        self.sys_name = path.split('/')[-1].split('.')[0]
        self.work_dir = path[0:-len(path.split('/')[-1])]

    @property
    def data_path(self):
        return self._data_path

    @data_path.setter
    def data_path(self, path):
        self._data_path = path
        self.sys_name = path.split('/')[-1].split('.')[1]
        self.work_dir = path[0:-len(path.split('/')[-1])]

    @property
    def com_range_list(self):
        return self._com_range_list

    @com_range_list.setter
    def com_range_list(self, a_list):
        self._com_range_list = [int(_) for _ in a_list]
        for (i, j) in zip(self.com_range_list[0::2], self.com_range_list[1::2]):
            self.com_range += list(range(i, j+1))

    def read_input(self):
        if not self.input_path:
            raise Exception('[ERROR] Input file not assigned!')
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'SPECIE':
                            self.specie_tags.append(line[1])
                            self.specie_numbers.append(int(line[2]))
                            self.total_species += 1
                        elif line[0] == 'COMRANGE':
                            self.com_range_list = [int(_) for _ in line[1:]]
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    # =================================================================================================================
    # PDB Methods.
    # =================================================================================================================
    def read_pdb(self):
        pdb_info = [' '.join(self.pdb_kw_list)]
        with open(self.pdb_path, 'r') as pdb_file:
            for line in pdb_file:
                if not line.startswith('ATOM'):
                    continue
                line = line.strip('\n')
                info = list()
                for kw in self.pdb_kw_list:
                    kw_range = self.pdb_kw_range[kw]
                    kw_info = line[kw_range[0]:kw_range[1]].strip(' ')
                    if kw_info:
                        info.append(kw_info)
                    else:
                        info.append('NULL')
                pdb_info.append(' '.join(info))
        self.pdb_df = pd.read_csv(StringIO('\n'.join(pdb_info)), delim_whitespace=True, dtype=self.pdb_kw_type)\
            .sort_values(by=['serial'])

    def data_to_pdb(self):
        new_pdb_path = '{}{}-new.pdb'.format(self.work_dir, self.sys_name)
        if not self.pdb_path:
            raise Exception('[ERROR] Reference .pdb file not given!'
                            'Info file from LAMMPS data file is not sufficient for generating a .pdb file!')
        self.read_pdb()
        self.read_data()
        pdb_part1 = pd.DataFrame.as_matrix(self.pdb_df, columns=['serial', 'name', 'altLoc', 'resName',
                                                                 'chainID', 'resSeq', 'iCode'])
        pdb_part2 = pd.DataFrame.as_matrix(self.data_atom_df, columns=['x', 'y', 'z'])
        pdb_part3 = pd.DataFrame.as_matrix(self.pdb_df, columns=['occupancy', 'tempFactor', 'element', 'charge'])
        new_pdb_info = np.concatenate((pdb_part1, pdb_part2, pdb_part3), axis=1)

        format_list = [6]
        for kw in self.pdb_kw_list:
            format_list.append(self.pdb_kw_range[kw][1] - self.pdb_kw_range[kw][0])
        with open(new_pdb_path, 'w') as out_file:
            out_file.write('HEADER\nTITLE\nREMARK\nREMARK\nREMARK\n')
            for line in new_pdb_info:
                line = [_ if _ != 'nan' else '' for _ in list(map(str, line))]

                out_line = 'ATOM{serial}  {name}{altLoc}{resName}{chainID}{resSeq}{iCode}{x}{y}{z}{occupancy}' \
                           '{tempFactor}{element}{charge}\n'\
                    .format(serial=line[0].rjust(7), name=line[1].ljust(3), altLoc=line[2].rjust(1),
                            resName=line[3].rjust(3), chainID=line[4].rjust(2),
                            resSeq=line[5].rjust(4), iCode=line[6].rjust(1),
                            x=line[7][:7].rjust(11), y=line[8][:7].rjust(8), z=line[9][:7].rjust(8),
                            occupancy=line[10].rjust(6), tempFactor=line[11].rjust(6),
                            element=line[12].rjust(12), charge=line[13].rjust(2))
                out_file.write(out_line)
            out_file.write('END\n')

    def pdb_find_com(self):
        self.read_input()
        self.read_pdb()
        com = np.array([0., 0., 0.])
        total_mass = 0.
        mass_dict = {}
        coord_list = []
        com_df = pd.DataFrame.as_matrix(self.pdb_df[self.pdb_df['serial'].isin(self.com_range)])
        for line in com_df:
            coord = line[7:10]
            coord_list.append(coord)
            element = line[12]
            if element not in mass_dict:
                mass_dict[element] = elem(element).mass
            atom_mass = mass_dict[element]
            total_mass += atom_mass
            com += np.array(atom_mass * coord, dtype=float)
        com /= total_mass
        distance_list = [np.linalg.norm(coord - com) for coord in coord_list]
        print('Center of Mass: {}'.format(com))
        print('Nearest Real Atom Coord: {}'.format(coord_list[int(np.argmin(distance_list))]))

    # =================================================================================================================
    # MOL2 Methods.
    # =================================================================================================================
    # Outdated function. Kept here for possible future use.
    # Use generate_mol2() instead of this.
    def generate_mol2_amber(self):
        trigger = input('System contains water/hydronium? [y/N]')
        if trigger == 'y':
            j_flag = '1'
        else:
            j_flag = '4'
        subprocess.call('$AMBERHOME/bin/antechamber -i %s -fi pdb -o %s.mol2 -fo mol2 -j %s -s 2 -rn' %
                        (self.pdb_path, self.sys_name, j_flag), shell=True)
        self.mol2_path = '{}{}.mol2'.format(self.work_dir, self.sys_name)

    def generate_mol2(self):
        mol2_coords = self.generate_coords()
        atom_tot = len(mol2_coords)
        mol2_bonds = self.generate_bonds()
        bond_tot = len(mol2_bonds)

        mol2_path = '{}{}.mol2'.format(self.work_dir, self.sys_name)
        with open(mol2_path, 'w') as mol2_file:
            mol2_file.write('@<TRIPOS>MOLECULE\nNAME\n')
            mol2_file.write('{} {} 1 0 0\n'.format(atom_tot, bond_tot))
            mol2_file.write('SMALL\nNO_CHARGES\n\n\n@<TRIPOS>ATOM\n')
            for line in mol2_coords:
                mol2_file.write('{:>7}  {:<3}{:>16.10}{:>11.10}{:>11.10} {:<3}{:>9}{:>4}{:>15}\n'
                                .format(line[0], line[1], line[2], line[3],
                                        line[4], line[5], line[6], line[7], line[8]))
            mol2_file.write('@<TRIPOS>BOND\n')
            for line in mol2_bonds:
                mol2_file.write(line)
            mol2_file.write('@<TRIPOS>SUBSTRUCTURE\n1 NAME 1 TEMP 0 **** **** 0 ROOT\n')
            self.mol2_path = '{}{}.mol2'.format(self.work_dir, self.sys_name)

    def generate_coords(self):
        coord_part1 = pd.DataFrame.as_matrix(self.pdb_df, columns=['serial', 'name', 'x', 'y', 'z'])
        coord_part3 = pd.DataFrame.as_matrix(self.pdb_df, columns=['resSeq', 'resName'])

        coord_part2 = list()
        coord_part4 = list()
        for number, specie in zip(self.specie_numbers, self.specie_tags):
            self.db.cursor.execute('SELECT MoleculeID FROM Molecules WHERE MoleculeType=\'{}\''.format(specie))
            molecule_id = self.db.cursor.fetchone()[0]
            self.db.cursor.execute('SELECT AtomType FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
            coord_part2.extend(self.db.cursor.fetchall() * number)
            self.db.cursor.execute('SELECT Charge FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
            coord_part4.extend(self.db.cursor.fetchall() * number)
        coord_part2 = np.array(coord_part2)
        coord_part4 = np.array(coord_part4)

        coord_info = np.concatenate((coord_part1, coord_part2, coord_part3, coord_part4), axis=1)
        return coord_info

    def generate_bonds(self):
        out_lines = list()
        bond_index = 1
        start = 0
        for i, specie in enumerate(self.specie_tags):
            self.db.cursor.execute('SELECT MoleculeID FROM Molecules WHERE MoleculeType=\'{}\''.format(specie))
            molecule_id = self.db.cursor.fetchone()[0]
            self.db.cursor.execute('SELECT Atom1ID, Atom2ID, BondType FROM Bonds WHERE MoleculeID={}'
                                   .format(molecule_id))
            bonds = self.db.cursor.fetchall()
            for number in range(self.specie_numbers[i]):
                for bond in bonds:
                    atom_1 = bond[0]
                    atom_2 = bond[1]
                    bond_type = bond[2]
                    out_lines.append('%s %s %s %s\n' % (bond_index, atom_1 + start, atom_2 + start, bond_type))
                    bond_index += 1
                self.db.cursor.execute('SELECT AtomID FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
                start += len(self.db.cursor.fetchall())
        return out_lines

    # =================================================================================================================
    # DATA Methods.
    # =================================================================================================================
    def read_data(self):
        atom_chunk = 'AtomID MolID AtomType Charge x y z\n'
        bond_chunk = 'BondID BondType Atom1ID Atom2ID\n'
        angle_chunk = 'AngleID AngleType Atom1ID Atom2ID Atom3ID\n'
        dihedral_chunk = 'DihedralID DihedralType Atom1ID Atom2ID Atom3ID Atom4ID\n'

        with open(self.data_path, 'r') as data_file:
            data_lines = np.array(data_file.readlines())
        atom_start = np.argwhere(data_lines == 'Atoms\n')[0][0]
        bond_start, bond = (np.argwhere(data_lines == 'Bonds\n')[0][0], True)\
            if np.argwhere(data_lines == 'Bonds\n') else (-1, False)
        angle_start, angle = (np.argwhere(data_lines == 'Angles\n')[0][0], True)\
            if np.argwhere(data_lines == 'Angles\n') else (-1, False)
        dihedral_start, dihedral = (np.argwhere(data_lines == 'Dihedrals\n')[0][0], True)\
            if np.argwhere(data_lines == 'Dihedrals\n') else (-1, False)
        atom_chunk += ''.join(data_lines[atom_start+1:bond_start])
        bond_chunk += ''.join(data_lines[bond_start+1:angle_start]) if bond else ''
        angle_chunk += ''.join(data_lines[angle_start+1:dihedral_start]) if angle else ''
        dihedral_chunk += ''.join(data_lines[dihedral_start+1:]) if dihedral else ''
        self.data_atom_df = pd.read_csv(StringIO(atom_chunk), delim_whitespace=True).sort_values(by=['AtomID'])
        self.data_bond_df = pd.read_csv(StringIO(bond_chunk), delim_whitespace=True).sort_values(by=['BondID'])
        self.data_angle_df = pd.read_csv(StringIO(angle_chunk), delim_whitespace=True).sort_values(by=['AngleID'])
        self.data_dihedral_df = pd.read_csv(StringIO(dihedral_chunk), delim_whitespace=True)\
            .sort_values(by=['DihedralID'])

    def modify_data(self):
        out_path = self.data_path
        in_path = out_path + '.ORIGINAL'
        subprocess.call('mv %s %s' % (out_path, in_path), shell=True)
        molecule_index = 1
        residue_ids = list()
        for i, specie in enumerate(self.specie_tags):
            for number in range(self.specie_numbers[i]):
                self.db.cursor.execute('SELECT MoleculeID FROM Molecules WHERE MoleculeType=\'{}\''.format(specie))
                molecule_id = self.db.cursor.fetchone()[0]
                self.db.cursor.execute('SELECT AtomID FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
                residue_ids.extend([molecule_index] * len(self.db.cursor.fetchall()))
                molecule_index += 1

        with open(in_path, 'r') as in_file:
            in_lines = in_file.readlines()
        skip_list = []
        with open(out_path, 'w') as out_file:
            for i, line in enumerate(in_lines):
                if i % 1000 == 0:
                    print('Done line %s' % i)
                if i in skip_list:
                    continue
                else:
                    lines = line.split()
                    if len(lines) == 0 or lines[0] != 'Atoms':
                        out_file.write(line)
                    else:
                        out_file.write(line)
                        j = 1
                        while True:
                            line = in_lines[i+j]
                            lines = line.split()
                            if len(lines) == 0:
                                out_file.write(line)
                                skip_list.append(i+j)
                                j += 1
                            elif lines[0] == 'Bonds':
                                break
                            else:
                                line = '%s %s %s %s %s %s %s\n' % (lines[0], residue_ids[int(lines[0])-1], lines[2],
                                                                   lines[3], lines[4], lines[5], lines[6])
                                out_file.write(line)
                                skip_list.append(i+j)
                                j += 1

    # =================================================================================================================
    # Miscellaneous Methods.
    # =================================================================================================================
    def generate_leaprc(self):
        leaprc_path = '%s.leaprc' % self.sys_name
        with open(leaprc_path, 'w') as leaprc_file:
            leaprc_file.write('source leaprc.gaff\n')
            leaprc_file.write('%s = loadmol2 %s\n' % (self.sys_name, self.mol2_path))
            leaprc_file.write('loadamberparams frcmod\n')
            leaprc_file.write('saveoff %s %s.lib\n' % (self.sys_name, self.sys_name))
            leaprc_file.write('saveamberparm %s %s.top %s.crd\n' % (self.sys_name, self.sys_name, self.sys_name))
            leaprc_file.write('quit')
        self.leaprc_path = leaprc_path

    def generate_top(self):
        in_path = self.pdb_path
        out_path = '%s.top' % self.sys_name
        with open(in_path, 'r') as in_file, open(out_path, 'w') as out_file:
            while True:
                line = in_file.readline()
                lines = line.split()
                if len(lines) == 0 or lines[0] != 'ATOM':
                    if lines[0] == 'END':
                        break
                    else:
                        continue
                else:
                    line_id = lines[1]
                    atom_name = lines[2]
                    molecule_type = lines[3]
                    self.db.cursor.execute('SELECT MoleculeID, KernelID FROM Molecules WHERE MoleculeType=\'{}\''
                                           .format(molecule_type))
                    molecule_id, kernel_id = self.db.cursor.fetchone()
                    self.db.cursor.execute('SELECT KernelTag FROM Atoms WHERE MoleculeID={} AND ATOMNAME=\'{}\''
                                           .format(molecule_id, atom_name))
                    kernel_tag = self.db.cursor.fetchone()[0]
                    if kernel_id:
                        to_write = ' '.join([line_id, str(kernel_id), str(kernel_tag), '\n'])
                    else:
                        to_write = ' '.join([line_id, '0', '0', '\n'])
                    out_file.write(to_write)

    def generate_sections(self):
        section_len = 0
        for i, molecule_type in enumerate(self.specie_tags):
            self.db.cursor.execute('SELECT MoleculeID FROM Molecules WHERE MoleculeType=\'{}\''.format(molecule_type))
            molecule_id = self.db.cursor.fetchone()[0]
            self.db.cursor.execute('SELECT AtomID FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
            section_len += len(self.db.cursor.fetchall()) * self.specie_numbers[i]
            self.sections[self.specie_tags[i]] = [self.specie_numbers[i], section_len]

