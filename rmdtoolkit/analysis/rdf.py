# Python 3.6.1

import numpy as np

from .base import AnalysisBase


class RDFAnalysis(AnalysisBase):
    def __init__(self, hist_range=(0., 20.), hist_bins=400, rho=None, **kwargs):
        super().__init__(coords=None, **kwargs)

        self.n_frame = 0
        self.n_center = None
        self.rho = rho

        self.hist_range = hist_range
        self.hist_bins = hist_bins
        self.bin_size = (hist_range[1] - hist_range[0]) / hist_bins

        self.bin_edges = None
        self.hist = np.zeros(hist_bins)

    def rdf(self, centers, particles):
        self.n_frame += 1

        if self.n_center is None:
            self.n_center = len(centers)
        else:
            if len(centers) != self.n_center:
                raise ValueError("{ClassName}.n_center is not consistent!".format(ClassName=self.__class__.__name__))

        if self.rho is None:
            self.rho = len(particles) / self.box_volume

        if self.bin_edges is None:
            self.bin_edges = np.histogram([0], bins=self.hist_bins, range=self.hist_range)[1]

        hists = np.array(list(map(lambda x: np.histogram(self.dist_pbc(x, particles, self.box_lengths),
                                                         bins=self.hist_bins, range=self.hist_range)[0],
                                  centers)))
        self.hist += np.sum(hists, axis=0)

    def normalize_rdf(self):
        self.hist /= 4 * np.pi * self.bin_edges[1:] ** 2 * self.bin_size * self.rho * self.n_center * self.n_frame

    @classmethod
    def dist_pbc(cls, origin, coords, box_lengths):
        dist = np.abs(coords - origin)
        pdist = box_lengths - dist
        mask = dist < box_lengths / 2
        dist = np.where(mask, dist, pdist)
        return np.linalg.norm(dist, axis=1)
