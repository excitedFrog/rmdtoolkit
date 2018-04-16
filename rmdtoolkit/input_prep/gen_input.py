# Python 3.6.1

import os
import subprocess

from rmdtoolkit.config.config_manager import ConfigManager
from rmdtoolkit.basic.model import Model
from rmdtoolkit.basic.result import Result

CM = ConfigManager()


class GenInput(Model, Result):
    def __init__(self):
        super().__init__()
        self.input_mol_path = str()

    def gen_data(self):
        if self.input_mol_path.endswith('psf'):
            self.psf_path = self.input_mol_path
            self.gen_data_charmm()
        elif self.input_mol_path.endswith('pdb'):
            self.pdb_path = self.input_mol_path
            self.gen_data_amber()
        elif self.input_mol_path.endswith('mol2'):
            self.mol2_path = self.input_mol_path
            self.gen_data_amber()
        else:
            raise Exception('[GenInput Error] Unrecognized input topology file extension.')

    def gen_data_charmm(self):  # TODO: implement this some time
        raise Exception('Feature not yet supported.')

    def gen_data_amber(self, rmd=False):
        tleap = CM.get('executables', 'tleap_executable')
        amb2lmp = CM.get('executables', 'amber2lammps')

        self.read_input()
        if self.input_mol_path.endswith('pdb'):
            print('Generating mol2 file......')
            self.read_pdb()
            self.generate_mol2()
        elif self.input_mol_path.endswith('mol2'):
            os.chdir(self.work_dir)

        print('Generating leaprc file......')
        self.generate_leaprc()
        print('Generating top file and crd file......')
        subprocess.call('%s -f %s' % (tleap, self.leaprc_path), shell=True)
        print('Generating data file......')
        subprocess.call('python2 %s' % amb2lmp, shell=True)
        print('Modifying data file......')
        self.data_path = 'data.%s' % self.sys_name
        self.modify_data()

        # Housekeeping
        subprocess.call('rm %s.top' % self.sys_name, shell=True)
        subprocess.call('rm leap.log', shell=True)
        subprocess.call('rm %s.crd' % self.sys_name, shell=True)
        subprocess.call('rm %s.leaprc' % self.sys_name, shell=True)
        subprocess.call('rm %s.lib' % self.sys_name, shell=True)

        if rmd:
            self.generate_top()
