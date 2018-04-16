# Python 3.6.1

import os

import numpy as np

from rmdtoolkit.basic.result import Result
from rmdtoolkit.basic.tool import pbc_dist


# wat_o_type and wat_h_type required, these attributes are located in Result() class
class WaterCluster(Result):
    def __init__(self):
        super().__init__()
        self.hyd_o_ids = list()
        self.hyd_h_ids = list()
        self.wat_o_ids = list()
        self.wat_h_ids = list()

        self.com_input = str()

        self.bin_size = 0.1
        self.num_bins = 5
        self.skip_first = 20

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'BIN_SIZE':
                            self.bin_size = float(line[1])
                        elif line[0] == 'BIN_NUM':
                            self.num_bins = int(line[1])
                        elif line[0] == 'SKIP_FIRST':
                            self.skip_first = int(line[1])
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def get_evb_clusters(self, max_shell=10):
        evb_clusters = list()
        center_coords = list()
        cluster_sizes = list()

        center_mol_ids = np.array(self.evb_info['ReactionCenters'])[:, 1]
        # TODO
        for center_mol_id in center_mol_ids:
            for line in self.trj_info:
                print(line)
                if int(line.split()[2]) == center_mol_id and int(line.split()[1]) in self.wat_o_type:
                    center_coords.append(np.array([float(i) for i in line.split()[3:6]]))

        for shell in self.evb_info['SHELLS']:
            cluster_sizes.append(len(shell))
            cluster_mol_ids = list()
            for line in shell:
                shell_id = line[2]
                if shell_id > max_shell:
                    break
                for mol_id in line[3:5]:
                    if mol_id not in cluster_mol_ids and mol_id != -1:
                        cluster_mol_ids.append(mol_id)
            evb_cluster = self.fetch_mol(cluster_mol_ids)
            evb_clusters.append(evb_cluster)

        evb_clusters = self.unfold_clusters(center_coords, evb_clusters)
        return evb_clusters, center_mol_ids, cluster_sizes

    def get_radius_clusters(self, radius=5.):
        radius_clusters = list()
        center_coords = list()
        center_mol_ids = list()
        for line in self.evb_info['ReactionCenters']:
            center_mol_ids.append(line[1])

        for center_mol_id in center_mol_ids:
            for line in self.trj_info:
                if int(line.split()[2]) == center_mol_id and int(line.split()[1]) in self.wat_o_type:
                    center_coords.append(np.array([float(i) for i in line.split()[3:6]]))

        box_len = [abs(bound[1] - bound[0]) for bound in self.bounds]
        for center_coord in center_coords:
            cluster_mol_ids = list()
            for line in self.trj_info:
                coord = np.array([float(i) for i in line.split()[3:6]])
                delta_r = pbc_dist(coord, center_coord, boxlen=box_len)
                if delta_r <= radius:
                    cluster_mol_id = int(line.split()[2])
                    if cluster_mol_id not in cluster_mol_ids:
                        cluster_mol_ids.append(cluster_mol_id)
            radius_cluster = self.fetch_mol(cluster_mol_ids)
            radius_clusters.append(radius_cluster)

        radius_clusters = self.unfold_clusters(center_coords, radius_clusters)
        return radius_clusters, center_mol_ids

    def unfold_clusters(self, center_coords, clusters):
        box_len = [abs(bound[1] - bound[0]) for bound in self.bounds]
        new_clusters = list()
        for center_coord, cluster in zip(center_coords, clusters):
            new_cluster = list()
            for atom in cluster:
                coord = np.array([float(i) for i in atom.split()[3:6]])
                new_coord = np.zeros(3)
                delta_r = coord - center_coord
                for i, ri in enumerate(delta_r):
                    if abs(ri) > box_len[i]/2:
                        new_coord[i] = coord[i] - np.sign(ri) * box_len[i]
                    else:
                        new_coord[i] = coord[i]
                new_atom = ' '.join([' '.join(atom.split()[0:3]), ' '.join(str(i) for i in new_coord)])
                new_cluster.append(new_atom)
            new_clusters.append(new_cluster)
        return new_clusters

    # proton sharing parameter, psp = |D_oa_h - D_od_h|
    def calculate_psp(self, cluster, center_mol_id):
        donor_o_coord = np.zeros(3)
        h_coords = list()
        acceptor_o_coords = list()

        for i, atom in enumerate(cluster):
            atom = atom.split()
            if int(atom[2]) == center_mol_id:
                if int(atom[1]) in self.wat_h_type:
                    h_coords.append(np.array([float(i) for i in atom[3:6]]))
                elif int(atom[1]) in self.wat_o_type:
                    donor_o_coord = np.array([float(i) for i in atom[3:6]])
            elif int(atom[2]) != center_mol_id and int(atom[1]) in self.wat_o_type:
                acceptor_o_coords.append(np.array([float(i) for i in atom[3:6]]))

        tmp_psps = list()
        for h_index, h_coord in enumerate(h_coords):
            od_h = np.linalg.norm(h_coord - donor_o_coord)
            tmp_oa_hs = list()
            for acceptor_o_coord in acceptor_o_coords:
                tmp_oa_h = np.linalg.norm(h_coord - acceptor_o_coord)
                tmp_oa_hs.append(tmp_oa_h)
            try:
                oa_h = min(tmp_oa_hs)
            except ValueError:
                oa_h = 100.0
                with open('Singleton', 'a') as f:
                    f.write('%s %s\n' % (center_mol_id, h_index))
            tmp_psp = abs(oa_h - od_h)
            tmp_psps.append(tmp_psp)
        psp = min(tmp_psps)
        return psp

    def psp_classify(self):
        os.chdir(self.work_dir)
        count = np.zeros(self.num_bins)
        self.find_file()
        self.open_file()

        for _ in range(self.skip_first):
            self.trj_read_frame()
        system_size = len(self.trj_info)
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            # Sometimes the last frame of trajectory is not written completely (due to 30 sec limit of slurm system)
            elif len(self.trj_info) != system_size:
                continue
            self.evb_match_trj()
            clusters, center_ids, cluster_sizes = self.get_evb_clusters()
            for cluster, center_id in zip(clusters, center_ids):
                delta = self.calculate_psp(cluster, center_id)
                i_bin = int(delta / self.bin_size)
                try:
                    count[i_bin] += 1
                except IndexError:
                    count[-1] += 1
        with open('{}{}.pspclass'.format(self.save_dir, self.save_tag), 'w') as save_file:
            for i, num in enumerate(count):
                save_file.write('{} {}\n'.format(i, num))
        self.close_file()

    def cation_classify(self):
        os.chdir(self.work_dir)
        self.find_file()
        self.open_file()

        cation_dict = dict()
        system_size = 0
        is_first_frame = True
        while True:
            self.trj_read_frame()
            if is_first_frame:
                system_size = len(self.trj_info)
            if self.trj_time == 'EOF':
                break
            # Sometimes the last frame of trajectory is not written completely (due to 30 sec limit of slurm system)
            elif len(self.trj_info) != system_size:
                continue
            self.evb_match_trj()
            num_water = len(self.evb_info['SHELLS'])
            if num_water not in cation_dict:
                cation_dict[num_water] = 0
            cation_dict[num_water] += 1
        with open('{}{}.cationclass'.format(self.save_dir, self.save_tag), 'w') as save_file:
            for num_water in cation_dict:
                save_file.write('{} {}\n'.format(num_water, cation_dict[num_water]))
        self.close_file()

    def save_cluster_xyz(self, cluster, save_tag):
        with open('%s%s.xyz' % (self.save_dir, save_tag), 'w') as f:
            f.write('%s\nNone\n' % len(cluster))
            for line in cluster:
                line = line.split()
                atom_type = line[1]
                atom_tag = self.atom_tag[atom_type]
                f.write('%s %s %s %s\n' % (atom_tag, line[3], line[4], line[5]))

    def log_struct_info(self, com_flag=False):
        with open('%s%s.dipole' % (self.save_dir, self.save_tag), 'a') as save_file:
            save_file.write('time r(H*-O*) r(Hw-Ow) r(Ow-O*) d(Ow-O*) r(COM-O*) d(COM-O*)\n')
            save_file.flush()
        os.chdir(self.work_dir)
        self.find_file()
        self.open_file()

        is_first_frame = True
        system_size = 0
        last_frame = 0
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            # Sometimes the last frame of trajectory is not written completely (due to 30 sec limit of slurm system)
            if not is_first_frame and len(self.trj_info) != system_size:
                continue
            if self.trj_time < last_frame:
                continue
            self.evb_match_trj()
            if not self.evb_info:
                break
            last_frame = self.trj_time

            if self.trj_time % 5000 == 0:
                print('{} READING FRAME {}'.format(self.save_tag, self.trj_time))

            clusters, center_mol_ids, cluster_sizes = self.get_evb_clusters(max_shell=1)
            com = np.zeros(3)
            if com_flag:
                com = self.find_com()
            for cluster, center_mol_id in zip(clusters, center_mol_ids):
                hydronium = np.zeros((4, 3))
                waters = list()
                water_mol_ids = list()
                hstar_count = 0
                hwat_counts = list()
                try:
                    for line in cluster:
                        line = line.split()
                        atom_type = int(line[1])
                        mol_id = int(line[2])
                        coordinate = [float(_) for _ in line[3:6]]
                        if mol_id == center_mol_id:  # Hydronium
                            if atom_type in self.wat_o_type:
                                hydronium[0][:] = coordinate
                            elif atom_type in self.wat_h_type:
                                hstar_count += 1
                                hydronium[hstar_count][:] = coordinate
                        else:
                            if mol_id not in water_mol_ids:
                                waters.append(np.zeros((3, 3)))
                                water_mol_ids.append(mol_id)
                                hwat_counts.append(0)
                            index = water_mol_ids.index(mol_id)
                            if atom_type in self.wat_o_type:
                                waters[index][0][:] = coordinate
                            elif atom_type in self.wat_h_type:
                                hwat_counts[index] += 1
                                # This is what is being "tried", sometimes RMD do not reassign mol_id
                                waters[index][hwat_counts[index]][:] = coordinate
                except IndexError:
                    continue
                r_com_ostar = np.zeros(3)
                d_com_ostar = 0.0
                if com_flag:
                    r_com_ostar = hydronium[0] - com
                    d_com_ostar = np.linalg.norm(r_com_ostar)
                r_hstars_ostar = hydronium[0] - hydronium[1:].sum(axis=0)/len(hydronium[1:])
                r_hwat_owat_list = list()
                r_owat_ostar_list = list()
                d_owat_ostar_list = list()
                for water in waters:
                    r_hwat_owat_list.append(water[0] - water[1:].sum(axis=0)/len(water[1:]))
                    r_owat_ostar = hydronium[0] - water[0]
                    r_owat_ostar_list.append(r_owat_ostar)
                    d_owat_ostar_list.append(np.linalg.norm(r_owat_ostar))
                with open('{}{}.dipole'.format(self.save_dir, self.save_tag), 'a') as save_file:
                    save_file.write('|'.join([str(self.trj_time),  # time
                                              str(center_mol_id),  # mol_id
                                              ' '.join([str(_) for _ in r_hstars_ostar]),  # r(H*-O*)
                                              ' '.join([str(_) for _ in [__ for __ in r_hwat_owat_list]]),  # r(Hw-Ow)
                                              ' '.join([str(_) for _ in [__ for __ in r_owat_ostar_list]]),  # r(Ow-O*)
                                              ' '.join([str(_) for _ in d_owat_ostar_list]),  # d(Ow-O*)
                                              ' '.join([str(_) if com_flag else com_flag for _ in r_com_ostar]),  # r(COM-O*)
                                              str(d_com_ostar) if com_flag else com_flag,
                                              '\n']))
                    save_file.flush()
        self.close_file()

    def calculate_vec_hstar_ostar(self):
        pass

    def log_struct_info_cluster(self, com_flag=False):
        with open('%s%s.dipole' % (self.save_dir, self.save_tag), 'a') as save_file:
            save_file.write('time r(H*-O*) r(Hw-Ow) r(Ow-O*) d(Ow-O*) r(COM-O*) d(COM-O*)\n')
            save_file.flush()
        os.chdir(self.work_dir)
        self.find_file()
        self.open_file()

        is_first_frame = True
        system_size = 0
        last_frame = 0
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            # Sometimes the last frame of trajectory is not written completely (due to 30 sec limit of slurm system)
            if not is_first_frame and len(self.trj_info) != system_size:
                continue
            if self.trj_time < last_frame:
                continue
            self.evb_match_trj()
            if not self.evb_info:
                break
            last_frame = self.trj_time

            if self.trj_time % 5000 == 0:
                print('{} READING FRAME {}'.format(self.save_tag, self.trj_time))
            clusters, center_mol_ids, cluster_sizes = self.get_evb_clusters(max_shell=1)
            com = np.zeros(3)
            if com_flag:
                com = self.find_com()
            for cluster, center_mol_id, cluster_size in zip(clusters, center_mol_ids, cluster_sizes):
                hydronium = np.zeros((4, 3))
                waters = list()
                water_mol_ids = list()
                hstar_count = 0
                hwat_counts = list()
                try:
                    for line in cluster:
                        line = line.split()
                        atom_type = int(line[1])
                        mol_id = int(line[2])
                        coordinate = [float(_) for _ in line[3:6]]
                        if mol_id == center_mol_id:  # Hydronium
                            if atom_type in self.wat_o_type:
                                hydronium[0][:] = coordinate
                            elif atom_type in self.wat_h_type:
                                hstar_count += 1
                                hydronium[hstar_count][:] = coordinate
                        else:
                            if mol_id not in water_mol_ids:
                                waters.append(np.zeros((3, 3)))
                                water_mol_ids.append(mol_id)
                                hwat_counts.append(0)
                            index = water_mol_ids.index(mol_id)
                            if atom_type in self.wat_o_type:
                                waters[index][0][:] = coordinate
                            elif atom_type in self.wat_h_type:
                                hwat_counts[index] += 1
                                # This is what is being "tried", sometimes RMD do not reassign mol_id
                                waters[index][hwat_counts[index]][:] = coordinate
                except IndexError:
                    continue
                r_com_ostar = np.zeros(3)
                d_com_ostar = 0.0
                if com_flag:
                    r_com_ostar = hydronium[0] - com
                    d_com_ostar = np.linalg.norm(r_com_ostar)
                r_hstars_ostar = hydronium[0] - hydronium[1:].sum(axis=0)/len(hydronium[1:])
                r_hwat_owat_list = list()
                r_owat_ostar_list = list()
                d_owat_ostar_list = list()
                for water in waters:
                    r_hwat_owat_list.append(water[0] - water[1:].sum(axis=0)/len(water[1:]))
                    r_owat_ostar = hydronium[0] - water[0]
                    r_owat_ostar_list.append(r_owat_ostar)
                    d_owat_ostar_list.append(np.linalg.norm(r_owat_ostar))
                with open('{}{}.dipole'.format(self.save_dir, self.save_tag), 'a') as save_file:
                    save_file.write('|'.join([str(self.trj_time),  # time
                                              str(center_mol_id),  # mol_id
                                              ' '.join([str(_) for _ in r_hstars_ostar]),  # r(H*-O*)
                                              ' '.join([str(_) for _ in [__ for __ in r_hwat_owat_list]]),  # r(Hw-Ow)
                                              ' '.join([str(_) for _ in [__ for __ in r_owat_ostar_list]]),  # r(Ow-O*)
                                              ' '.join([str(_) for _ in d_owat_ostar_list]),  # d(Ow-O*)
                                              ' '.join([str(_) if com_flag else com_flag for _ in r_com_ostar]),  # r(COM-O*)
                                              str(d_com_ostar) if com_flag else com_flag,
                                              str(cluster_size),
                                              '\n']))
                    save_file.flush()
        self.close_file()

    def log_struct_info_necessary(self, com_flag=False):
        with open('%s%s.dipole' % (self.save_dir, self.save_tag), 'a') as save_file:
            save_file.write('time r(H*-O*) r(Hw-Ow) r(Ow-O*) d(Ow-O*) r(COM-O*) d(COM-O*)\n')
            save_file.flush()
        os.chdir(self.work_dir)
        self.find_file()
        self.open_file()

        is_first_frame = True
        system_size = 0
        last_frame = 0
        while True:
            self.trj_read_frame()
            if self.trj_time == 'EOF':
                break
            # Sometimes the last frame of trajectory is not written completely (due to 30 sec limit of slurm system)
            if not is_first_frame and len(self.trj_info) != system_size:
                continue
            if self.trj_time < last_frame:
                continue
            self.evb_match_trj()
            if not self.evb_info:
                break
            last_frame = self.trj_time

            if self.trj_time % 5000 == 0:
                print('{} READING FRAME {}'.format(self.save_tag, self.trj_time))

            clusters, center_mol_ids, cluster_sizes = self.get_evb_clusters(max_shell=1)
            com = np.zeros(3)
            if com_flag:
                com = self.find_com()
            for cluster, center_mol_id in zip(clusters, center_mol_ids):
                hydronium = np.zeros((4, 3))
                hstar_count = 0
                for line in cluster:
                    line = line.split()
                    atom_type = int(line[1])
                    mol_id = int(line[2])
                    coordinate = [float(_) for _ in line[3:6]]
                    if mol_id == center_mol_id:  # Hydronium
                        if atom_type in self.wat_o_type:
                            hydronium[0][:] = coordinate
                        elif atom_type in self.wat_h_type:
                            hstar_count += 1
                            hydronium[hstar_count][:] = coordinate

                r_com_ostar = np.zeros(3)
                d_com_ostar = 0.0
                if com_flag:
                    r_com_ostar = hydronium[0] - com
                    d_com_ostar = np.linalg.norm(r_com_ostar)
                r_hstars_ostar = hydronium[0] - hydronium[1:].sum(axis=0)/len(hydronium[1:])
                with open('{}{}.dipole'.format(self.save_dir, self.save_tag), 'a') as save_file:
                    save_file.write('|'.join([str(self.trj_time),  # time
                                              str(center_mol_id),  # mol_id
                                              ' '.join([str(_) for _ in r_hstars_ostar]),  # r(H*-O*)
                                              ' '.join([str(_) if com_flag else com_flag for _ in r_com_ostar]),  # r(COM-O*)
                                              str(d_com_ostar) if com_flag else com_flag,
                                              '\n']))
                    save_file.flush()
        self.close_file()
