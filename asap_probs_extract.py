import numpy as np
import numba

class distanced_extractor():
    def __init__(self, src):
        self.hist = dict()
        self.source = src

    def build(self, thresh):