# Python 3.6.1

import sys
import subprocess
import os
import random
import time
import itertools

import imaplib

from config.config_manager import ConfigManager
from basic.tool import text_inform

CM = ConfigManager()

GEN_FILES_DIR = CM.get('basics', 'gen_files_dir')
LAMMPS_EXE = CM.get('lammps', 'lammps_executable')
# LAMMPS_INPUT_DIR stores STATIC files that is required for every simulation. e.g. something.par
LAMMPS_INPUT_DIR = CM.get('lammps', 'lammps_input_dir')

QOS = CM.get('slurm', 'special_partition_prefix')
PARTITION_LIST = CM.get_list('slurm', 'partition_list')
DEFAULT_PARTITION = CM.get('slurm', 'default_partition')

EMAIL_SERVER = CM.get('email', 'email_server')
EMAIL_USER_NAME = CM.get('email', 'email_account_name')
EMAIL_PWD = CM.get('email', 'email_account_pwd')
EMAIL_ADDRESS = CM.get('email', 'email_address')

MY_PHONE_NUMBER = CM.get('twilio', 'phone_number')


class JobMonitor(object):
    def __init__(self):
        super().__init__()
        self.monitor = dict()
        self.idle_partition = str()
        self.iter_lists = list()
        self.prefix_dir = str()

    def run(self, is_restart):
        if not is_restart:
            for iter_elements in itertools.product(*self.iter_lists):
                system_tag = iter_elements[0]  # may not be [0], but it gotta be somewhere in iter_elements right?
                self.lammps_submit(iter_elements, system_tag, is_init=True)
        while True:
            time.sleep(900)
            self.update_monitor()
            sys.stdout.write('Job status monitor update successful {}\n'.format(time.ctime()))
            sys.stdout.flush()
            for job_name in self.monitor:
                if self.monitor[job_name][1] == 'TIMEOUT':
                    iter_elements = job_name.split('-')
                    system_tag = iter_elements[0]
                    self.lammps_submit(iter_elements, system_tag, is_init=False)
                elif self.monitor[job_name][1] == 'FAILED':
                    text_message = 'Job %s has failed.' % job_name
                    text_inform(MY_PHONE_NUMBER, text_message)
                    self.monitor[job_name][1] = 'FAILED-INFORMED'

    def find_idle_partition(self):
        idle_partition = str()
        idle_amount = 0
        for partition in PARTITION_LIST:
            try:
                info = subprocess.check_output('sinfo -p %s|grep idle' % partition, shell=True)
            except subprocess.CalledProcessError:
                continue
            info = info.split()
            if int(info[3].decode()) > idle_amount:
                idle_amount = int(info[3].decode())
                idle_partition = info[0].decode()
        if idle_partition == str():
            idle_partition = DEFAULT_PARTITION
        self.idle_partition = idle_partition.strip('*')

    def lammps_submit(self, iter_elements, system_tag, is_init):
        """
        This function is to-some-degree ad hoc. Modify before use.
        :param iter_elements: tuple or list of dir names, no '/' sign; e.g. ($TEMP, $SYSTEM_SIZE, $FORCE_FIELD, etc.)
        :param system_tag: useful for identifying system name.
        :param is_init: True for initial run(minimization, assign vel, etc.), False for restarting.
        self.prefix_dir/
            |-systemA/            <- "data_dir/"
            |      |-in.data
            |      |   |-Trial1/    <- "work_dir/"
            |      |   |-Trial2/    <- "work_dir/"
            |      |   |-Trial3/    <- "work_dir/"
        ...     ...     ...
        lmp_input/                <- $LAMMPS_INPUT_DIR
            |-ForceField.par
            |-EVB.pra
        ...     ...     ...
        """
        # Make working directories.
        data_dir = '%s%s/' % (self.prefix_dir, system_tag)  # stores in.data files for different simulations
        work_dir = '%s%s/' % (self.prefix_dir, '/'.join(iter_elements))  # where the simulation takes place
        subprocess.call('mkdir -p %s' % work_dir, shell=True)

        # Find previous runs, determine output basename
        day_count = 1
        file_list = os.listdir(work_dir)
        for file_path in file_list:
            if file_path.endswith('lammpstrj'):
                day_count += 1

        # Generate submit.sh
        job_name = '-'.join(item for item in iter_elements)
        self.find_idle_partition()
        with open('%ssubmit.sh' % work_dir, 'w') as submit:
            submit.write('#!/usr/bin/env sh\n\n')
            if self.idle_partition.startswith(QOS):
                submit.write('#SBATCH --qos=%s\n' % QOS)
                submit.write('#SBATCH --time=24:00:00\n')
            else:
                submit.write('#SBATCH --time=1-12:00:00\n')
            submit.write('#SBATCH --partition=%s\n' % self.idle_partition)
            submit.write('#SBATCH --job-name=%s\n' % job_name)
            submit.write('#SBATCH --output=%s.out\n' % day_count)
            submit.write('#SBATCH --mail-user=%s\n' % EMAIL_ADDRESS)
            with open('%slammps_submit_gen' % GEN_FILES_DIR, 'r') as f:
                for line in f:
                    if line.startswith('exe='):
                        line.strip('\n')
                        line += LAMMPS_EXE
                        line += '\n'
                    submit.write(line)

        # Generate in.lmp.sh
        with open('%sin.lmp.sh' % work_dir, 'w') as inlmp:
            inlmp.write('variable  count  string  %s\n' % day_count)

            if is_init:
                inlmp.write('variable  seed  string %s\n' % random.randint(1, 10000000))

            with open('%s%s.lmp_gen' % (GEN_FILES_DIR, system_tag), 'r') as f:
                for line in f:
                    inlmp.write(line)
                    line = line.split()
                    if len(line) > 0 and line[0] == 'neigh_modify':
                        if is_init:
                            inlmp.write('read_data  in.data\n')
                        else:
                            inlmp.write('read_restart  1.restart\n')

                    elif len(line) > 0 and line[0] == 'thermo_style':
                        if is_init:
                            inlmp.write('minimize 1.0e-4 1.0e-6 1000 10000\n')
                            inlmp.write('reset_timestep  0\n')
                            inlmp.write('velocity  all create ${Temp_simu} ${seed} rot yes dist gaussian\n')

        # Copy other essential files to working directory.
        if is_init:
            subprocess.call('cp %s/* %s' % (LAMMPS_INPUT_DIR, work_dir), shell=True)
            subprocess.call('cp %s/* %s' % (data_dir, work_dir), shell=True)

        info = subprocess.check_output('sbatch submit.sh', shell=True)
        info = info.decode().strip('\n').split()
        job_id = info[-1]
        os.chdir(self.prefix_dir)

        sys.stdout.write('{} submitted as {}'.format(job_name, job_id))
        self.monitor[job_name] = [job_id, 'in_queue']

    def update_monitor(self):
        status_list = ['Began', 'FAILED', 'TIMEOUT', 'CANCELLED', 'COMPLETED']

        connection = imaplib.IMAP4_SSL(EMAIL_SERVER, '993')
        connection.login(EMAIL_USER_NAME, EMAIL_PWD)
        connection.select(readonly=0)
        retcode, unread_mails = connection.search(None, '(UNSEEN)')

        for mail_id in unread_mails[0].split():
            typ, data = connection.fetch(mail_id, '(RFC822)')
            content = data[0][1].decode()

            job_id = 'Null'
            job_name = 'Null'
            job_status = 'Null'

            for job_status in status_list:
                flg = content.find(job_status)
                if flg != -1:
                    break

            content = content.split('\r\n')
            for i in range(len(content)):
                line = content[i]
                if line.startswith('Subject:'):
                    line = line.split()
                    job_id = line[2][7:]
                    job_name = line[3][5:]
            self.monitor[job_name] = [job_id, job_status]
            connection.uid('STORE', mail_id, '+FLAGS', '/SEEN')

        connection.close()
        connection.logout()
