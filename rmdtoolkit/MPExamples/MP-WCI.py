# Python 3.6.1

import multiprocessing as mp

from rmdtoolkit.result_analysis.interface import Interface


def analysis_worker(input_path):
    cls = Interface()
    cls.input_path = input_path
    cls.wci_worker()


with mp.Pool as Pool:
    Pool.starmap(None, None)
