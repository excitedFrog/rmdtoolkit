# Python 3.6.1

import sqlite3 as sql
import pandas as pd
from io import StringIO

from mendeleev import element

from rmdtoolkit.config.config_manager import ConfigManager

CM = ConfigManager()
DATABASE_PATH = CM.get('basics', 'database_path')


class MyDB(object):
    def __init__(self):
        super().__init__()
        self.db_path = DATABASE_PATH
        self.mol_file = str()
        self.mol_df = pd.DataFrame()
        self.connection = None
        self.cursor = None
        self.initialize()
        self.pdb_kw_list = ['serial', 'name', 'altLoc', 'resName', 'chainID', 'resSeq', 'iCode', 'x', 'y', 'z',
                            'occupancy', 'tempFactor', 'element', 'charge']
        self.pdb_kw_range = {'serial': [6, 11], 'name': [13, 16], 'altLoc': [16, 17], 'resName': [17, 20],
                             'chainID': [21, 22], 'resSeq': [22, 26], 'iCode': [26, 27],
                             'x': [30, 38], 'y': [38, 46], 'z': [46, 54],
                             'occupancy': [54, 60], 'tempFactor': [60, 66], 'element': [76, 78], 'charge': [78, 80]}
        self.pdb_kw_type = {'serial': int, 'name': str, 'altLoc': str, 'resName': str, 'chainID': str, 'resSeq': int,
                            'iCode': str, 'x': float, 'y': float, 'z': float, 'occupancy': str, 'tempFactor': str,
                            'element': str, 'charge': str}

    def initialize(self):
        self.connection = sql.connect(self.db_path)
        self.cursor = self.connection.cursor()

    def drop_all(self):
        self.cursor.execute('DROP TABLE Molecules')
        self.cursor.execute('DROP TABLE Atoms')
        self.cursor.execute('DROP TABLE Elements')
        self.cursor.execute('DROP TABLE Bonds')
        self.connection.commit()

    def construct(self):
        self.cursor.execute('CREATE TABLE Molecules (\n'
                            '  MoleculeID INTEGER NOT NULL PRIMARY KEY,\n'
                            '  MoleculeType VARCHAR(255),\n'
                            '  MoleculeDesc VARCHAR(255),\n'
                            '  KernelID INTEGER\n'
                            ')')
        self.cursor.execute('CREATE TABLE Atoms (\n'
                            '  AtomID INTEGER NOT NULL,\n'
                            '  MoleculeID INTEGER NOT NULL,\n'
                            '  ElementID INTEGER NOT NULL,\n'
                            '  AtomName VARCHAR(255) NOT NULL,\n'
                            '  AtomType VARCHAR(255) NOT NULL,\n'
                            '  Charge REAL NOT NULL,\n'
                            '  KernelTag INTEGER,\n'
                            '  PRIMARY KEY (MoleculeID, AtomID),\n'
                            '  FOREIGN KEY (MoleculeID) REFERENCES Molecules(MoleculeID),\n'
                            '  FOREIGN KEY (ElementID) REFERENCES Elements(ElementID)\n'
                            ')')
        self.cursor.execute('CREATE TABLE Elements (\n'
                            '  ElementID INTEGER NOT NULL PRIMARY KEY,\n'
                            '  ElementName VARCHAR(255) NOT NULL,\n'
                            '  ElementMass REAL NOT NULL\n'
                            ')')
        self.cursor.execute('CREATE TABLE Bonds (\n'
                            '  BondID INTEGER NOT NULL,\n'
                            '  MoleculeID INTEGER NOT NULL,\n'
                            '  BondType VARCHAR(255),\n'
                            '  Atom1ID INTEGER NOT NULL,\n'
                            '  Atom2ID INTEGER NOT NULL,\n'
                            '  PRIMARY KEY (MoleculeID, BondID),\n'
                            '  FOREIGN KEY (MoleculeID) REFERENCES Molecules(MoleculeID)\n'
                            ')')
        self.connection.commit()

    def gen_element_table(self):
        for i in range(1, 119):
            symbol = element(i).symbol
            mass = element(i).mass
            self.cursor.execute('INSERT INTO Elements(ElementID, ElementName, ElementMass) VALUES ({}, \'{}\', {})'
                                .format(i, symbol, mass))
            self.connection.commit()

    # Function for migrating info from a dict from an older version of this package
    def dict_to_db(self, dictionary):
        self.cursor.execute('SELECT * FROM Molecules')
        molecule_id = len(self.cursor.fetchall())
        for molecule_type in dictionary:
            molecule_id += 1
            self.cursor.execute('INSERT INTO Molecules(MoleculeID, MoleculeType, MoleculeDesc)\n'
                                'VALUES ({}, \'{}\', \'{}\')'.format(molecule_id, molecule_type, ' '))
            self.connection.commit()

            bond_id = 0
            for atom1id in dictionary[molecule_type]['Bonds']:
                for atom2id in dictionary[molecule_type]['Bonds'][atom1id]:
                    bond_type = dictionary[molecule_type]['Bonds'][atom1id][atom2id]
                    bond_id += 1
                    self.cursor.execute('INSERT INTO Bonds(BondID, MoleculeID, BondType, Atom1ID, Atom2ID)\n'
                                        'VALUES ({}, {}, \'{}\', {}, {})'
                                        .format(bond_id, molecule_id, bond_type, atom1id, atom2id))
                    self.connection.commit()

            for atom_id in dictionary[molecule_type]['Atoms']:
                info = dictionary[molecule_type]['Atoms'][atom_id]
                element_name = info[1]
                self.cursor.execute('SELECT ElementID FROM Elements WHERE ElementName=\'{}\''.format(element_name))
                element_id = self.cursor.fetchone()[0]
                atom_name = info[0]
                atom_type = info[2]
                charge = info[3]
                self.cursor.execute('INSERT INTO Atoms(AtomID, MoleculeID, ElementID, AtomName, AtomType, Charge)\n'
                                    'VALUES ({}, {}, {}, \'{}\', \'{}\', {})'
                                    .format(atom_id, molecule_id, element_id, atom_name, atom_type, charge))
                self.connection.commit()

    def read_pdb(self):
        pdb_info = [' '.join(self.pdb_kw_list)]
        with open(self.mol_file, 'r') as pdb_file:
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
        self.mol_df = pd.read_csv(StringIO('\n'.join(pdb_info)), delim_whitespace=True, dtype=self.pdb_kw_type)\
            .sort_values(by=['serial'])

    def import_info_pdb(self, molecule_id, molecule_type, molecule_desc='', kernel_id=0):
        self.cursor.execute('SELECT MoleculeID FROM Molecules')
        molecule_ids = [_[0] for _ in self.cursor.fetchall()]
        self.cursor.execute('SELECT MoleculeType FROM Molecules')
        molecule_types = [_[0] for _ in self.cursor.fetchall()]
        if molecule_id in molecule_ids or molecule_type in molecule_types:
            raise Exception('[DATABASE ERROR] Given MoleculeID/MoleculeType already exists!')

        self.cursor.execute('INSERT INTO Molecules(MoleculeID, MoleculeType, MoleculeDesc, KernelID)\n'
                            'VALUES ({}, \'{}\', \'{}\', {})'
                            .format(molecule_id, molecule_type, molecule_desc, kernel_id))
        self.read_pdb()
        atom_info = pd.DataFrame.as_matrix(self.mol_df, columns=['serial', 'name', 'element'])
        for row in atom_info:
            atom_id = row[0]
            self.cursor.execute('SELECT ElementID FROM Elements WHERE ElementName=\'{}\''.format(row[2]))
            element_id = self.cursor.fetchone()[0]
            atom_name = row[1]
            atom_type = ''
            charge = 0.
            self.cursor.execute('INSERT INTO Atoms(AtomID, MoleculeID, ElementID, AtomName, AtomType, Charge)\n'
                                'VALUES ({}, {}, {}, \'{}\', \'{}\', {})'
                                .format(atom_id, molecule_id, element_id, atom_name, atom_type, charge))
        self.connection.commit()

    def delete_molecule(self, molecule_id=None, molecule_type=None):
        if not molecule_id and not molecule_type:
            raise Exception('[ERROR] MoleculeID or MoleculeType for the molecule to be deleted is needed.')
        if not molecule_id:
            self.cursor.execute('SELECT MoleculeID FROM Molecules WHERE MoleculeType=\'{}\''.format(molecule_type))
            molecule_id = self.cursor.fetchone()[0]
        self.cursor.execute('DELETE FROM Molecules WHERE MoleculeID={}'.format(molecule_id))
        self.cursor.execute('DELETE FROM Atoms WHERE MoleculeID={}'.format(molecule_id))
        self.cursor.execute('DELETE FROM Bonds WHERE MoleculeID={}'.format(molecule_id))
        self.connection.commit()


if __name__ == '__main__':
    db = MyDB()
    db.mol_file = '/Users/zheful/Projects/rmdtoolkit/water.pdb'
    # db.delete_molecule(molecule_type='WAT')
    db.import_info_pdb(1, 'WAT', molecule_desc='Water, H2O')
