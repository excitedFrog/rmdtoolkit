# Python 3.6.1

import numpy as np
from collections import defaultdict


class Graph(object):
    def __init__(self):
        self.graph = defaultdict(list)

    def add_edge(self, u, v):
        self.graph[u].append(v)

    def add_edges(self, edges):
        for edge in edges:
            self.add_edge(edge[0], edge[1])
            self.add_edge(edge[1], edge[0])

    def bfs(self, s, shell_lim=np.inf):
        visited = [False] * len(self.graph)
        queue = list()
        queue.append(s)
        visited[s] = True
        i_shell = 0

        shells = list()
        nodes = list()
        while queue:
            s = queue.pop(0)
            shells.append(i_shell)
            nodes.append(s)
            if not queue:
                i_shell += 1
                if i_shell > shell_lim:
                    break
            for i in self.graph[s]:
                if not visited[i]:
                    queue.append(i)
                    visited[i] = True
        return nodes, shells
