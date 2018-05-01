# Python 3.6.1

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from rmdtoolkit.basic.result import Result


class Interface(Result):
    def __init__(self):
        super().__init__()
        self.module_handle = 'INTERFACE'

        self.wci_bounds = np.zeros((3, 2))
        self.num_grid = None
        self.mesh_grid = np.array([])
        self.bandwidth = 2.4
        self.den_values = None
        self.xyz = {'x': 0, 'y': 1, 'z': 2}
        self.interface = None

        self.d_interface = 0.016
        self.interval = 0.5
        self._mode = 'complete'
        self.supported_modes = ['complete', 'fast']

        self.den_atoms = list()  # Density of which is used to determine the interface
        self.pmf_atoms = list()  # PMF atom
        self.dist_results = list()

        #Temp
        self.zcec = float()

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, string):
        if string in self.supported_modes:
            self._mode = string
        else:
            raise Exception('[{}] Unrecognized calculation mode \'{}\'. Aborting.' .format(self.module_handle, string))

    def read_input(self):
        super().read_input()
        with open(self.input_path, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if len(line) > 0:
                    try:
                        if line[0] == 'WCI_X':
                            self.wci_bounds[0] = np.array(list(map(float, line[1:3])))
                        elif line[0] == 'WCI_Y':
                            self.wci_bounds[1] = np.array(list(map(float, line[1:3])))
                        elif line[0] == 'WCI_Z':
                            self.wci_bounds[2] = np.array(list(map(float, line[1:3])))
                        elif line[0] == 'INTERFACE_DENSITY':
                            self.d_interface = float(line[1])
                        elif line[0] == 'MESH_INTERVAL':
                            self.interval = float(line[1])
                        elif line[0] == 'MODE':
                            self.mode = line[1]
                        elif line[0] == 'DEN_ATOM':
                            self.den_atoms.append(line[1:])
                        elif line[0] == 'PMF_ATOM':
                            self.pmf_atoms.append(line[1:])
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def gaussian_kde(self, origin):
        diff = np.abs(self.den_values - origin)  # if no pbc
        pdiff = self.box_len - diff  # if all pbc
        mask = diff < self.box_len/2  # find which element need pbc
        diff = np.where(mask, diff, pdiff)  # mask'em
        tdiff = np.dot(diff, 1 / self.bandwidth)
        energy = np.sum(tdiff * tdiff, axis=1) / 2.
        return np.sum(np.exp(-energy)) / (2 * np.pi * self.bandwidth**2)

    def gaussian_kde_new(self):  # Not in use, since somehow this fully vectorized method is slower than map. Funny.
        diffs = np.abs(self.den_values[None, ...] - self.mesh_grid[:, None, :])
        pdiffs = self.box_len - diffs
        masks = diffs < self.box_len / 2
        diffs = np.where(masks, diffs, pdiffs)
        tdiffs = np.dot(diffs, 1 / self.bandwidth)
        energies = np.sum(tdiffs * tdiffs, axis=2) / 2.
        return np.sum(np.exp(-energies), axis=1) / (2 * np.pi * self.bandwidth**2)

    def dist_to_interface(self, origin):
        diff = np.abs(self.interface - origin)
        pdiff = self.box_len - diff
        mask = diff < self.box_len/2
        diff = np.where(mask, diff, pdiff)
        dist = np.min(np.linalg.norm(diff, axis=1))
        return dist

    def z_to_interface(self, origin):
        return np.abs(
            origin[2] - self.interface[np.argmin(np.linalg.norm(self.interface[:, 0:2] - origin[0:2], axis=1)), 2])

    def wci(self):  # Willard-Chandler Interface
        if not self.mesh_grid.any():
            self.wci_bounds = np.where(self.wci_bounds.astype(bool), self.wci_bounds, self.bounds)
            self.num_grid = (np.array(list(map(lambda _: _[1]-_[0], self.wci_bounds))) / self.interval).astype(int)
            x, y, z = np.mgrid[self.wci_bounds[0][0]:self.wci_bounds[0][1]:self.num_grid[0]*1j,
                               self.wci_bounds[1][0]:self.wci_bounds[1][1]:self.num_grid[1]*1j,
                               self.wci_bounds[2][0]:self.wci_bounds[2][1]:self.num_grid[2]*1j]
            self.mesh_grid = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

        self.den_values = self.atoms_find(self.den_atoms)
        density_grid = np.abs(np.array(list(map(self.gaussian_kde, self.mesh_grid))) - self.d_interface)\
            .reshape(self.num_grid)
        mesh_grid = self.mesh_grid.reshape(np.append(self.num_grid, 3))

        if self.mode == 'fast':
            self.interface = mesh_grid[np.where(density_grid < 0.002)].T
        elif self.mode == 'complete':
            index = np.argmin(density_grid, axis=2)
            xs, ys = np.mgrid[0:self.num_grid[0], 0:self.num_grid[1]]
            self.interface = mesh_grid[xs, ys, index].reshape((self.num_grid[0]*self.num_grid[1], 3))

    def save_wci(self):
        with open('{}{}.wci'.format(self.save_dir, self.save_tag), 'a') as save_file:
                save_file.write('# WCI at {}\n'.format(self.trj_time))
                np.savetxt(save_file, self.interface, fmt='%.4f')

    def wci_pmf(self):
        self.wci()
        origins = self.atoms_find(self.pmf_atoms)
        self.dist_results = list(map(self.dist_to_interface, origins))
        # Temp
        self.zcec = self.find_cec()[0][2]

    def wci_zpmf(self):
        self.wci()
        origins = self.atoms_find(self.pmf_atoms)
        self.dist_results = list(map(self.z_to_interface, origins))
        # Temp
        self.zcec = self.find_cec()[0][2]

    # This does not calculate PMF directly, but logs the COLVAR(dist to interface) instead.
    # The PMF shall be calculated with WHAM, as there is a series of biased trajectories.
    def save_wci_pmf(self):
        with open('{}{}.wci-pmf'.format(self.save_dir, self.save_tag), 'a') as save_file:
            for dist in self.dist_results:
                save_file.write('1.000 {} {}\n'.format(dist, self.zcec))

    def wci_worker(self):
        self.analysis_template(inner_compute_func=self.wci,
                               inner_save_func=self.save_wci,
                               outer_compute_func=self.void_func,
                               outer_save_func=self.void_func)

    def wci_plot_worker(self):
        self.analysis_template(inner_compute_func=self.wci,
                               inner_save_func=self.plot,
                               outer_compute_func=self.void_func,
                               outer_save_func=self.void_func)

    def wci_pmf_worker(self):
        self.analysis_template(inner_compute_func=self.wci_pmf,
                               inner_save_func=self.save_wci_pmf,
                               outer_compute_func=self.void_func,
                               outer_save_func=self.void_func)

    def wci_zpmf_worker(self):
        self.analysis_template(inner_compute_func=self.wci_zpmf,
                               inner_save_func=self.save_wci_pmf,
                               outer_compute_func=self.void_func,
                               outer_save_func=self.void_func)

    def plot(self):
        print(self.frame_tot)
        self.interface = self.interface.T
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if self._mode == 'fast':
            ax.scatter(self.interface[0], self.interface[1], self.interface[2])
        if self._mode == 'complete':
            xs = self.interface[0].reshape((self.num_grid[0], self.num_grid[1])).T[0]
            ys = self.interface[1].reshape((self.num_grid[0], self.num_grid[1]))[0]
            xs, ys = np.meshgrid(xs, ys)
            zs = self.interface[2].reshape((self.num_grid[0], self.num_grid[1]))

            ax.plot_surface(xs, ys, zs, cmap='jet')
        ax.auto_scale_xyz([-12, 12], [-12, 12], [-30, -15])
        ax.set_xlabel('X coord')
        ax.set_ylabel('Y coord')
        ax.set_zlabel('Z coord')
        plt.savefig('frame{}.jpg'.format(self.frame_tot))
        plt.close()
        self.interface = self.interface.T
