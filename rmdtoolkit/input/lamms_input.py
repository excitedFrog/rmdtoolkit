# Python 3.6.1

import os

from .base import InputBase
from ..parser import PDBParser


class DataInput(InputBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.pdb_path = os.path.join(self.work_dir, '{}.pdb'.format(self.base_name))
        self.mol2_path = os.path.join(self.work_dir, '{}.mol2'.format(self.base_name))
        self.frcmod_path = os.path.join(self.work_dir, '{}.frcmod'.format(self.base_name))
        self.leaprc_path = os.path.join(self.work_dir, '{}.leaprc'.format(self.base_name))
        self.data_path = os.path.join(self.work_dir, '{}.data'.format(self.base_name))

    def gen_data_amber(self):
        tleap = self.cm.get('executables', 'tleap_executable')
        amb2lmp = self.cm.get('executables', 'amber2lammps')

    def gen_leaprc(self):
        with open(self.leaprc_path, 'w') as leaprc_file:
            leaprc_file.write('source leaprc.gaff\n'
                              '{base_name} = loadmol2 {mol2_path}\n'
                              'loadamberparams {frcmod_path}\n'
                              'saveoff {base_name} {base_name}.lib\n'
                              'saveamberparm %s %s.top %s.crd\n'
                              'quit'
                              .format(base_name=self.base_name, mol2_path=self.mol2_path, frcmod_path=self.frcmod_path))

