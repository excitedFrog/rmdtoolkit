import numpy as np

from .base import AnalysisBase


class CoordNumAnalysis(AnalysisBase):
    def __init__(self):
        super().__init__(coords=None)
        self.center_coords = np.array([])
        self.dist_coords = np.array([])

    def coord_number_worker(self, center, cutoff, mode):
        dist_list = self.dist_pbc(center, self.dist_coords, self.box_lengths)
        if mode == 'hard':
            coord_num = np.sum(dist_list < cutoff)
        elif mode == 'poly':
            coord_num = np.sum(list(map(lambda dist: self.polynomial_coord_num(dist, cutoff), dist_list)))
        else:
            raise ValueError("{ClassName} Unrecognized mode!".format(ClassName=self.__class__.__name__))
        return coord_num

    def coord_number(self, cutoff=2.5, mode='poly'):
        return np.array(list(map(lambda center: self.coord_number_worker(center, cutoff, mode), self.center_coords)))

    @classmethod
    def dist_pbc(cls, origin, coords, box_lengths):
        dist = np.abs(coords - origin)
        pdist = box_lengths - dist
        mask = dist < box_lengths / 2
        dist = np.where(mask, dist, pdist)
        return np.linalg.norm(dist, axis=1)

    # White and Voth, JCTC, 2014 10 3023-3030.
    @classmethod
    def polynomial_coord_num(cls, d, cutoff, pow1=6, pow2=12, w=0.3):
        return 1 if d <= cutoff else (1 - ((d - cutoff) / w) ** pow1) / (1 - ((d - cutoff) / w) ** pow2)
