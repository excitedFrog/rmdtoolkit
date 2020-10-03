# Python 3.6.1

import numpy as np

from .base import AnalysisBase


class Assign(AnalysisBase):
    def __init__(self, o_coords, h_coords, **kwargs):
        super().__init__(o_coords, **kwargs)
        self.o_coords = o_coords
        self.h_coords = h_coords

    # def assign_oh(self):

