# Copyright (C) Aleksei Ponomarenko-Timofeev
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import numpy as np
import pairdata as pd
from scipy.constants import speed_of_light
from sklearn.metrics import r2_score
from auxfun import *

__author__ = 'Aleksei Ponomarenko-Timofeev'

def check_chan_NF(chan: pd.chan = None, freq: float = 60e9):
    src_coords = chan.src.coords
    dst_coords = chan.dest.coords

    for i in chan.paths.items():
        if i[1].interactions.__len__() > 0:
            if (np.linalg.norm(i[1].interactions[0].coords - src_coords) < (10.0 * speed_of_light / freq)) or \
                (np.linalg.norm(i[1].interactions[-1].coords - dst_coords) < (10.0 * speed_of_light / freq)):
                i[1].near_field_failed = True
            else:
                k = False
                for j in i[1].interactions:
                    if not k:
                        prev_node = j
                        k = True
                        pass
                    elif np.linalg.norm(np.asarray(prev_node.coords) - np.asarray(j.coords)) < (10.0 * speed_of_light / freq):
                        i[1].near_field_failed = True
                        break

    fixed_pow = 0.0

    for i in chan.paths.items():
        if i[1].near_field_failed:
            pass
        else:
            fixed_pow += i[1].pow

    chan.valid_pow = fixed_pow

def check_data_NF(data: pd.data_stor = None, freq: float = 60e9):
    for i in data.txs.items():
        for j in i[1].chans_to_pairs.items():
            check_chan_NF(j[1], freq=freq)


class cluster(list):
    def __init__(self, threshold: float = 0.2, inipath: pd.path = None):
        assert inipath is not None, "Cluster needs an initial path to build around!"
        list.__init__(self)
        self.inipath = inipath
        self.threshold = threshold
        self.pow = 0.0
        self.intcoords = inipath.chan.src.coords
        for i in inipath.interactions:
            [self.intcoords.append(j) for j in i.coords]
        [self.intcoords.append(j) for j in inipath.chan.dest.coords]

    def append(self, path: pd.path):
        assert path is not None, "Can't append None to cluster! Check Your code!"
        self.pow += path.pow
        list.append(self, path)

    def cappend(self, path : pd.path):
        assert path is not None, "Can't append None to cluster! Check Your code!"

        pathcoords = path.chan.src.coords
        for i in path.interactions:
            [pathcoords.append(j) for j in i.coords]
        [pathcoords.append(j) for j in path.chan.dest.coords]

        if 1.0 - r2_score(self.intcoords, pathcoords) > self.threshold:
            return False
        else:
            self.pow += path.pow
            list.append(self, path)
            return True