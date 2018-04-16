# Python 3.6.1

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from rmdtoolkit.basic.result import Result
from rmdtoolkit.basic.tool import pbc_dist
from scipy import stats


class Interface(Result):
    def __init__(self):
        super().__init__()
        self.module_handle = 'INTERFACE'

        self.wci_bounds = np.zeros((3, 2))
        self.num_grid = None
        self.bandwidth = 2.4
        self.values = None
        self.xyz = {'x': 0, 'y': 1, 'z': 2}
        self.interface = None

        self.d_interface = 0.016
        self.interval = 0.5
        self._mode = 'complete'
        self.mode_supported = ['complete', 'fast']

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, string):
        if string in self.mode_supported:
            self.mode = string
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
                    except IndexError:
                        raise Exception('No value found after tag \'{}\' in input file.'.format(line[0]))

    def gaussian_kde(self, origin):
        diff = np.abs(self.values - origin)  # if no pbc
        pdiff = self.box_len - diff  # if all pbc
        mask = diff < self.box_len/2  # find which element need pbc
        diff = np.where(mask, diff, pdiff)  # mask'em
        tdiff = np.dot(diff, 1 / self.bandwidth)
        energy = np.sum(tdiff * tdiff, axis=1) / 2.
        return np.sum(np.exp(-energy)) / (2 * np.pi * self.bandwidth**2)

    def willard_chandler_interface(self):
        self.wci_bounds = np.where(self.wci_bounds.astype(bool), self.wci_bounds, self.bounds)
        self.num_grid = (np.array(list(map(lambda _: _[1]-_[0], self.wci_bounds))) / self.interval).astype(int)
        x, y, z = np.mgrid[self.wci_bounds[0][0]:self.wci_bounds[0][1]:self.num_grid[0]*1j,
                           self.wci_bounds[1][0]:self.wci_bounds[1][1]:self.num_grid[1]*1j,
                           self.wci_bounds[2][0]:self.wci_bounds[2][1]:self.num_grid[2]*1j]
        mesh_grid = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

        self.values = self.atoms_find()
        density_grid = np.abs(np.array(list(map(self.gaussian_kde, mesh_grid))) - self.d_interface).reshape(self.num_grid)
        mesh_grid = mesh_grid.reshape(np.append(self.num_grid, 3))

        if self.mode == 'fast':
            self.interface = mesh_grid[np.where(density_grid < 0.002)].T
        elif self.mode == 'complete':
            index = np.argmin(density_grid, axis=2)
            xs, ys = np.mgrid[0:self.num_grid[0], 0:self.num_grid[1]]
            self.interface = mesh_grid[xs, ys, index].reshape((self.num_grid[0]*self.num_grid[1], 3)).T

    def worker(self):
        self.analysis_worker_template(self.willard_chandler_interface, self.none_func, self.save_wci)

    def save_wci(self):
        pass

    def plot(self, mode='complete'):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if mode == 'fast':
            ax.scatter(self.interface[0], self.interface[1], self.interface[2])
        if mode == 'complete':
            xs = self.interface[0].reshape((self.num_grid[0], self.num_grid[1])).T[0]
            ys = self.interface[1].reshape((self.num_grid[0], self.num_grid[1]))[0]
            xs, ys = np.meshgrid(xs, ys)
            zs = self.interface[2].reshape((self.num_grid[0], self.num_grid[1]))

            ax.plot_surface(xs, ys, zs, cmap='jet')
        ax.auto_scale_xyz([-12, 12], [-12, 12], [-30, -15])
        ax.set_xlabel('X coord')
        ax.set_ylabel('Y coord')
        ax.set_zlabel('Z coord')
        plt.savefig('grid3d.jpg')


