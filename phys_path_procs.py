# Copyright (C) Aleksei Ponomarenko-Timofeev
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__author__ = 'Aleksei Ponomarenko-Timofeev'

import numpy as np
from pairdata import DataStorage, Channel, Path
from scipy.constants import speed_of_light
from sklearn.metrics import r2_score
from auxfun import *


def check_chan_NF(chan: Channel = None, freq: float = 60e9):
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
                    elif np.linalg.norm(np.asarray(prev_node.coords) -
                                        np.asarray(j.coords)) < (10.0 * speed_of_light / freq):
                        i[1].near_field_failed = True
                        break

    fixed_pow = 0.0

    for i in chan.paths.items():
        if i[1].near_field_failed:
            pass
        else:
            fixed_pow += i[1].pow

    chan.valid_pow = fixed_pow


def check_data_nf(data: DataStorage = None, freq: float = 60e9):
    for i in data.txs.items():
        for j in i[1].chans_to_pairs.items():
            check_chan_NF(j[1], freq=freq)


class Cluster(list):
    def __init__(self, threshold: float = 1e-5, inipath: Path = None):
        assert inipath is not None, 'Cluster needs an initial path to build around!'
        list.__init__(self)
        self.inipath = inipath
        self.threshold = threshold
        self.pow = inipath.pow
        self.inivect = inipath.chan.src.coords.tolist()
        self.inipath.chan.clusters[inipath] = self
        self.inipath.cluster = self
        inipath.cluster = self
        for i in inipath.interactions:
            for j in i.coords.tolist():
                self.inivect.append(j)

        for j in inipath.chan.dest.coords.tolist():
            self.inivect.append(j)
        list.append(self, inipath)

    def append(self, path: Path):
        assert path is not None, 'Can\'t append None to cluster! Check Your code!'
        self.pow += path.pow
        path.cluster = self
        list.append(self, path)
        self.inipath.cluster = self
        # Recalculate

    def cappend(self, path: Path):
        assert path is not None, 'Can\'t append None to cluster! Check Your code!'
        if path.cluster is None:
            pathcoords = path.chan.src.coords.tolist()
            for i in path.interactions:
                for j in i.coords.tolist():
                    pathcoords.append(j)
            for j in path.chan.dest.coords.tolist():
                pathcoords.append(j)

            if not (pathcoords.__len__() == self.inivect.__len__()):
                return False

            if (1.0 - r2_score(self.inivect, pathcoords)) > self.threshold:
                return False

            path.cluster = self
            self.pow += path.pow
            self.inipath.cluster = self
            list.append(self, path)
            return True

        return False


def gen_chan_clusters(chan: Channel = None, threshold: float = 1e-5, nff: bool = False):
    assert chan is not None, 'None given instead of channel in gen_clusters!'
    for i in chan.paths.items():
        if not nff:
            if i[1].cluster is None:
                c = Cluster(threshold=threshold, inipath=i[1])
                for j in chan.paths.items():
                    if nff and not j[1].near_field_failed:
                        c.cappend(j[1])
                    elif not nff:
                        c.cappend(j[1])
        elif not i[1].near_field_failed and (i[1].cluster is None):
            c = Cluster(threshold=threshold, inipath=i[1])
            for j in chan.paths.items():
                if nff and not j[1].near_field_failed:
                    c.cappend(j[1])
                elif not nff:
                    c.cappend(j[1])


def gen_data_clusters(data: DataStorage, threshold: float = 1e-5, nff: bool = False):
    for i in data.txs.items():
        for j in i[1].chans_to_pairs.items():
            gen_chan_clusters(j[1], threshold=threshold, nff=nff)
