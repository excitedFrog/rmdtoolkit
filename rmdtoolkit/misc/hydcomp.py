# Python 3.6.1

import numpy as np
import itertools
from collections import defaultdict


class HydComp(object):
    def __init__(self, directed=False):
        self.box_len = float()
        self.graph = defaultdict(set)
        self.h_alloc = np.array([])
        self.o_coords = np.array([])
        self.h_coords = np.array([])
        self.directed = directed

        self.shell = -1
        self.shell_list = np.array([])
        self.node_list = np.array([])
        self.visited_nodes = set()

    def add_connection(self, node1, node2):
        self.graph[node1].add(node2)
        if not self.directed:
            self.graph[node2].add(node1)

    def add_connections(self, connections):
        for node1, node2 in connections:
            self.add_connection(node1, node2)

    def is_connected(self, node1, node2):
        return node1 in self.graph and node2 in self.graph[node1]

    def __str__(self):
        return '{}({})'.format(self.__class__.__name__, dict(self.graph))

    def dfs_util(self, node, visited):
        self.shell += 1
        visited[node] = True
        self.visited_nodes.add(node)
        self.node_list = np.append(self.node_list, node)
        self.shell_list = np.append(self.shell_list, [self.shell])

        to_visit = self.graph[node] - self.visited_nodes
        if not to_visit:
            self.shell = 0
            return 0
        for i in to_visit:
            if not visited[i]:
                self.dfs_util(i, visited)

    def search_shell_dfs(self, start_node):
        self.shell = -1
        self.shell_list = list()
        self.node_list = list()
        self.visited_nodes = set()

        visited = dict.fromkeys(list(self.graph.keys()), False)
        self.dfs_util(start_node, visited)
        shells = self.shell_list
        nodes = self.node_list

        return shells, nodes

    def d_oo(self, max_shell=1):
        shells = np.array([])
        d_oos = np.array([])
        for i in np.arange(max_shell):
            shell = i + 1
            inner_list = self.node_list[np.where(self.shell_list == i)]
            outer_list = self.node_list[np.where(self.shell_list == i + 1)]
            for inner in inner_list:
                deltas = list()
                doos = list()
                for outer in outer_list:
                    if outer not in self.graph[inner]:
                        continue

                    o_coord1 = self.o_coords[int(inner)]
                    o_coord2 = self.o_coords[int(outer)]
                    h_coords = self.h_coords[np.append(self.h_alloc[int(inner)], self.h_alloc[int(outer)])]

                    ds = np.abs(np.array(list(map(lambda _: self.dist_pbc(origin=o_coord1, coord=_,
                                                                          box_lengths=self.box_len), h_coords)))
                                - np.array(list(map(lambda _: self.dist_pbc(origin=o_coord2, coord=_,
                                                                            box_lengths=self.box_len), h_coords))))
                    doo = self.dist_pbc(o_coord1, o_coord2, self.box_len)
                    doos.append(doo)
                    deltas.append(min(ds))
                if not deltas:
                    continue
                doos = np.array(doos)
                doo = doos[np.argmin(deltas)]
                d_oos = np.append(d_oos, [doo])
                shells = np.append(shells, [shell])
        return shells, d_oos

    def angle_oho(self, max_shell=1):  # only valid when max_shell=1
        shells = np.array([])
        angles = np.array([])
        for i in np.arange(max_shell):
            shell = i + 1
            inner_list = self.node_list[np.where(self.shell_list == i)]
            outer_list = self.node_list[np.where(self.shell_list == i + 1)]
            for inner in inner_list:
                deltas = list()
                angs = list()
                for outer in outer_list:
                    if outer not in self.graph[inner]:
                        continue
                    o_coord1 = self.o_coords[int(inner)]
                    o_coord2 = self.o_coords[int(outer)]
                    h_coords = self.h_coords[np.append(self.h_alloc[int(inner)], self.h_alloc[int(outer)])]

                    ds = np.abs(np.array(list(map(lambda _: self.dist_pbc(origin=o_coord1, coord=_,
                                                                          box_lengths=self.box_len), h_coords)))
                                - np.array(list(map(lambda _: self.dist_pbc(origin=o_coord2, coord=_,
                                                                            box_lengths=self.box_len), h_coords))))
                    shared_h_coord = h_coords[np.argmin(ds)]
                    angle = self.angle_between(shared_h_coord, o_coord1, o_coord2, box_lengths=self.box_len)

                    deltas.append(min(ds))
                    angs.append(angle)
                if not deltas:
                    continue
                angs = np.array(angs)
                angle = angs[np.argmin(deltas)]
                angles = np.append(angles, [angle])
                shells = np.append(shells, [shell])
        return shells, angles

    def calculate_sp_struct(self):  # sp is for "special pair"
        hyd = self.node_list[np.where(self.shell_list == 0)]
        hyd_o_coord = self.o_coords[int(hyd)]
        wats = self.node_list[np.where(self.shell_list == 1)]
        wat_o_coords = self.o_coords[wats.astype(int)]

        if wat_o_coords.size == 0:
            return np.nan, np.nan, np.nan

        dist_oo = np.array(list(map(lambda _: self.dist_pbc(origin=hyd_o_coord, coord=_, box_lengths=self.box_len),
                                    wat_o_coords)))
        arg = np.argmin(dist_oo)
        doo = dist_oo[arg]
        wat = wats[arg]
        wat_o_coord = wat_o_coords[arg]

        h_coords = self.h_coords[np.append(self.h_alloc[int(hyd)], self.h_alloc[int(wat)])]
        dist_oh = np.array(list(map(lambda _: self.dist_pbc(origin=hyd_o_coord, coord=_, box_lengths=self.box_len),
                                    h_coords))) \
            + np.array(list(map(lambda _: self.dist_pbc(origin=wat_o_coord, coord=_, box_lengths=self.box_len),
                                h_coords)))
        shared_h_coord = h_coords[np.argmin(dist_oh)]

        angle = self.angle_between(shared_h_coord, wat_o_coord, hyd_o_coord, box_lengths=self.box_len)

        delta = np.abs(self.dist_pbc(origin=shared_h_coord, coord=hyd_o_coord, box_lengths=self.box_len) -
                       self.dist_pbc(origin=shared_h_coord, coord=wat_o_coord, box_lengths=self.box_len))
        return delta, doo, angle

    def calculate_delta(self):
        deltas = np.array([])
        hyd = self.node_list[np.where(self.shell_list == 0)]
        hyd_o_coord = self.o_coords[int(hyd)]
        wats = self.node_list[np.where(self.shell_list == 1)]
        if wats.size == 0:
            return 10.

        for wat in wats:
            wat_o_coord = self.o_coords[int(wat)]
            h_coords = self.h_coords[np.append(self.h_alloc[int(hyd)], self.h_alloc[int(wat)])]
            dist = np.array(list(map(lambda _: self.dist_pbc(origin=hyd_o_coord, coord=_, box_lengths=self.box_len),
                                     h_coords))) \
                + np.array(list(map(lambda _: self.dist_pbc(origin=wat_o_coord, coord=_, box_lengths=self.box_len),
                                    h_coords)))
            shared_h_coord = h_coords[np.argmin(dist)]

            delta = np.abs(self.dist_pbc(origin=shared_h_coord, coord=hyd_o_coord, box_lengths=self.box_len) -
                           self.dist_pbc(origin=shared_h_coord, coord=wat_o_coord, box_lengths=self.box_len))
            deltas = np.append(deltas, delta)
        return np.min(deltas)

    def calculate_hyd_dipole(self):
        hyd = self.node_list[np.where(self.shell_list == 0)]
        o_coord = self.o_coords[int(hyd)]
        h_coords = self.h_coords[self.h_alloc[int(hyd)]]
        h_coords = np.array(list(map(lambda coord: self.vec_pbc(o_coord, coord, self.box_len), h_coords)))
        dipole = np.array([0, 0, 0]) - np.sum(h_coords, axis=0) / 3
        return dipole

    def calculate_theta(self, ref):
        dipole = self.calculate_hyd_dipole()
        dipole /= np.linalg.norm(dipole)
        ref /= np.linalg.norm(ref)
        angle = np.arccos(np.clip(np.dot(dipole, ref), -1.0, 1.0))
        return angle

    def return_type(self):
        o_count = len(self.node_list)
        h_count = len(np.concatenate(self.h_alloc[self.node_list.astype(int)]))
        return "O%sH%s" % (o_count, h_count)

    @classmethod
    def dist_pbc(cls, origin, coord, box_lengths):
        dist = np.abs(coord - origin)
        pdist = box_lengths - dist
        mask = dist < box_lengths / 2
        dist = np.where(mask, dist, pdist)
        return np.linalg.norm(dist)

    @classmethod
    def vec_pbc(cls, origin, coord, box_lengths):
        vec = coord - origin
        pvec = vec - box_lengths * np.sign(vec)
        mask = abs(vec) < box_lengths / 2
        return np.where(mask, vec, pvec)

    @classmethod
    def angle_between(cls, c, p1, p2, box_lengths):
        v1 = cls.vec_pbc(c, p1, box_lengths)
        v2 = cls.vec_pbc(c, p2, box_lengths)
        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
