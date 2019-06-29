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

import sqlite3
import mysql.connector as msqlc
from mysql.connector.constants import ClientFlag
from numpy.linalg import norm
from numpy import asarray
from numpy import zeros
import concurrent.futures as cof
from multiprocessing import cpu_count
import scipy.io as sio
from sql_adapters import *
from siso_sql import *
from auxfun import l2db


def _load_paths_txthread(self, txsids, host, user, pw, dbname):
    dbconn = msqlc.connect(host=host, user=user, password=pw, client_flags=[ClientFlag.SSL], database=dbname)

    dbcurs = dbconn.cursor()
    dbcurs.execute(TX_PAIRST.format(min(txsids), max(txsids), min(txsids), max(txsids)))
    txp = dbcurs.fetchall()

    for i in txp:
        dst = self.rxs[i[-2]]
        self.txs[i[-1]].chans_to_pairs[dst] = Channel(dest=dst, src=self.txs[i[-1]])
        self.rxs[i[-2]].chans_to_pairs[self.txs[i]] = self.txs[i[-1]].chans_to_pairs[dst]
        self.txs[i[-1]].chans_to_pairs[dst].pow = i[1] * 1e3
        self.txs[i[-1]].chans_to_pairs[dst].delay = i[2]
        self.txs[i[-1]].chans_to_pairs[dst].delay_spread = i[3]
        self.txs[i[-1]].chans_to_pairs[dst].dist = norm(self.txs[i].coords - self.rxs[i[-2]].coords)
        self.txs[i[-1]].chans_to_pairs[dst].chid = i[0]

    for i in txsids:
        for j in self.txs[i].chans_to_pairs.keys():
            dbcurs.execute(CHAN_PTH.format(self.txs[i].chans_to_pairs[j].chid))
            d = dbcurs.fetchall()
            d = sorted(d, key=lambda t: t[1], reverse=True)
            for k in d[0:(self.npaths if d.__len__() >= self.npaths else d.__len__())]:
                self.txs[i].chans_to_pairs[j].paths[k[0]] = Path()
                self.txs[i].chans_to_pairs[j].paths[k[0]].chan = self.txs[i].chans_to_pairs[j]
                self.txs[i].chans_to_pairs[j].paths[k[0]].pathid = k[0]
                self.txs[i].chans_to_pairs[j].paths[k[0]].pow = k[1] * 1e3
                self.txs[i].chans_to_pairs[j].paths[k[0]].fspl = k[7]
                self.txs[i].chans_to_pairs[j].paths[k[0]].phase = k[8]
                self.txs[i].chans_to_pairs[j].paths[k[0]].delay = k[2]
                self.txs[i].chans_to_pairs[j].paths[k[0]].aod = k[3]
                self.txs[i].chans_to_pairs[j].paths[k[0]].eod = k[4]
                self.txs[i].chans_to_pairs[j].paths[k[0]].aoa = k[5]
                self.txs[i].chans_to_pairs[j].paths[k[0]].eoa = k[6]
    dbconn.close()


def _load_paths_rxthread(self, rxsids, host, user, pw, dbname):
    dbconn = msqlc.connect(host=host, user=user, password=pw, client_flags=[ClientFlag.SSL], database=dbname)

    dbcurs = dbconn.cursor()
    dbcurs.execute(RX_PAIRST.format(min(rxsids), max(rxsids), min(rxsids), max(rxsids)))
    rxp = dbcurs.fetchall()

    for i in rxp:
        dst = self.rxs[i[-1]]
        dst.chans_to_pairs[self.txs[i[-2]]] = Channel(src=self.txs[i[-2]], dest=dst)
        self.txs[i[-2]].chans_to_pairs[dst] = dst.chans_to_pairs[self.txs[i[-2]]]
        dst.chans_to_pairs[self.txs[i[-2]]].pow = i[1] * 1e3
        dst.chans_to_pairs[self.txs[i[-2]]].delay = i[2]
        dst.chans_to_pairs[self.txs[i[-2]]].delay_spread = i[3]
        dst.chans_to_pairs[self.txs[i[-2]]].dist = norm(dst.coords - self.txs[i[-2]].coords)
        dst.chans_to_pairs[self.txs[i[-2]]].chid = i[0]

    for i in rxsids:
        for j in self.rxs[i].chans_to_pairs.keys():
            dbcurs.execute(CHAN_PTH.format(self.rxs[i].chans_to_pairs[j].chid))
            d = dbcurs.fetchall()
            d = sorted(d, key=lambda t: t[1], reverse=True)
            for k in d[0:(self.npaths if d.__len__() >= self.npaths else d.__len__())]:
                self.rxs[i].chans_to_pairs[j].paths[k[0]] = Path()
                self.rxs[i].chans_to_pairs[j].paths[k[0]].chan = self.rxs[i].chans_to_pairs[j]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].pathid = k[0]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].pow = k[1] * 1e3
                self.rxs[i].chans_to_pairs[j].paths[k[0]].fspl = k[7]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].phase = k[8]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].delay = k[2]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].aod = k[3]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].eod = k[4]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].aoa = k[5]
                self.rxs[i].chans_to_pairs[j].paths[k[0]].eoa = k[6]

    dbconn.close()


def _load_iters_txs(self, txsids, host, user, pw, dbname, store):
    dbconn = msqlc.connect(host=host, user=user, password=pw,
                           client_flags=[ClientFlag.SSL], database=dbname)
    dbcurs = dbconn.cursor()

    for i in txsids:
        for j in self.txs[i].chans_to_pairs.items():
            dbcurs.execute(INTERS_SPEC_CHAN.format(j[1].chid))
            inters = dbcurs.fetchall()

            for k in j[1].paths.items():
                dists = list()
                precoords = j[1].src.coords
                for l in inters:
                    if l[4] == k[1].pathid:
                        if store:
                            intr = Interaction()
                            intr.path = k[1]
                            intr.coords = asarray([l[0], l[1], l[2]])
                            intr.typ = l[3]
                            dists.append(norm(precoords - intr.coords))
                            precoords = intr.coords
                        else:
                            intr = False
                        k[1].interactions.append(intr)
                dists.append(norm(precoords - j[1].dest.coords))
                k[1].length = sum(dists)

    dbconn.close()


def _load_iters_rxs(self, rxsids, host, user, pw, dbname, store):
    dbconn = msqlc.connect(host=host, user=user, password=pw,
                           client_flags=[ClientFlag.SSL], database=dbname)
    dbcurs = dbconn.cursor()

    for i in rxsids:
        for j in self.rxs[i].chans_to_pairs.items():
            dbcurs.execute(INTERS_SPEC_CHAN.format(j[1].chid))
            inters = dbcurs.fetchall()
            for k in j[1].paths.items():
                dists = list()
                precoords = j[1].src.coords
                for l in inters:
                    if l[4] == k[1].pathid:
                        if store:
                            intr = Interaction()
                            intr.path = k[1]
                            intr.coords = asarray([l[0], l[1], l[2]])
                            intr.typ = l[3]
                            dists.append(norm(precoords - intr.coords))
                            precoords = intr.coords
                        else:
                            intr = False
                        k[1].interactions.append(intr)
                dists.append(norm(precoords - j[1].dest.coords))
                k[1].length = sum(dists)

    dbconn.close()


class DataStorage:
    def __init__(self, conf: str = None, threaded: bool = True, dbname: str = ''):
        self.dbname = dbname
        self.pw = None
        self.uname = None
        self.dbconn = None
        self.dbcurs = None
        self.possible_inters = dict()
        self.threshold_chan = -115
        self.threshold_path = -125
        self.threaded = threaded
        self.npaths = 250

        if conf is not None:
            conff = open(conf)
            self.host = conff.readline().strip('\n')
            self.user = conff.readline().strip('\n')
            self.pw = conff.readline().strip('\n')

            print('Connecting to {} as {}'.format(self.host, self.user))
            conff.close()
            self.dbconn = msqlc.connect(host=self.host, user=self.user, password=self.pw,
                                        client_flags=[ClientFlag.SSL], database=self.dbname)
            self.dbcurs = self.dbconn.cursor()
            self.dbcurs.execute('USE {};'.format(self.dbname))
            self.dbname = self.dbname
        else:
            self.dbconn = sqlite3.connect(self.dbname)
            self.dbcurs = self.dbconn.cursor()

        self.txs = TXConnector(dbcursor=self.dbcurs, master=self)
        self.rxs = RXConnector(dbcursor=self.dbcurs, master=self)


    def __repr__(self):
        return 'SQL data storage'

    def __str__(self):
        if hasattr(self, 'host'):
            return 'Database {} at {}.'.format(self.dbname, self.host)
        else:
            return 'Database in file {}.'.format(self.dbname)

    def __link(self):
        if self.dbconn is None:
            if not hasattr(self, 'host'):
                self.dbconn = sqlite3.connect(self.dbname)
                self.dbcurs = self.dbconn.cursor()
            else:
                self.dbconn = msqlc.connect(host=self.host, user=self.user, password=self.pw,
                                            client_flags=[ClientFlag.SSL], database=self.dbname)
                self.dbcurs = self.dbconn.cursor()
                self.dbcurs.execute('USE {};'.format(self.dbname))
                self.dbname = self.dbname

    def __unlink(self):
        if self.dbconn is not None:
            self.dbconn.close()
            self.dbcurs = None
            self.dbconn = None

    def load_rxtx(self):
        print('Loading TX/RX nodes...', end='', flush=True)

        self.__link()
        # Read TX data
        self.dbcurs.execute(TX_EXTR)
        j = self.dbcurs.fetchall()

        for i in j:
            n = Node('TX')
            n.coords = asarray([i[1], i[2], i[3]])
            n.node_id = i[0]
            n.setid = i[4]
            self.txs[i[0]] = n

        # Read RX data
        self.dbcurs.execute(RX_EXTR)
        j = self.dbcurs.fetchall()

        for i in j:
            n = Node('RX')
            n.coords = asarray([i[1], i[2], i[3]])
            n.node_id = i[0]
            n.setid = i[4]
            self.rxs[i[0]] = n

        print('Success!', flush=True)

        self.__unlink()

    def load_paths(self, npaths: int = 250):
        print('Loading paths...', end='', flush=True)
        self.npaths = npaths

        self.__link()

        if hasattr(self, 'host') and self.threaded:
            txs_per_thread = self.txs.__len__() / self.nthreads
            rxs_per_thread = self.rxs.__len__() / self.nthreads

            thread_pool = cof.ThreadPoolExecutor(max_workers=self.nthreads)

            if txs_per_thread > rxs_per_thread:
                print('TXward...', end='', flush=True)
                txs_per_thread = int(txs_per_thread)
                txs = list()

                for i in self.txs.keys():
                    txs.append(i)
                    if txs.__len__() == txs_per_thread:
                        thread_pool.submit(_load_paths_txthread, self, txs, self.host, self.user, self.pw,
                                           self.dbname)
                        txs = list()

                if txs.__len__() > 0 or txs_per_thread <= 1:
                    thread_pool.submit(_load_paths_txthread, self, txs, self.host, self.user, self.pw, self.dbname)
            else:
                print('RXward...', end='', flush=True)
                rxs_per_thread = int(rxs_per_thread)
                rxs = list()

                for i in self.rxs.keys():
                    rxs.append(i)
                    if rxs.__len__() == rxs_per_thread:
                        thread_pool.submit(_load_paths_rxthread, self, rxs, self.host, self.user, self.pw,
                                           self.dbname)
                        rxs = list()

                if rxs.__len__() > 0 or rxs_per_thread <= 1:
                    thread_pool.submit(_load_paths_rxthread, self, rxs, self.host, self.user, self.pw, self.dbname)

            thread_pool.shutdown()
        else:
            for i in self.txs.keys():
                self.dbcurs.execute(TX_PAIRS.format(i, i))
                k = self.dbcurs.fetchall()
                for j in k:
                    dst = self.rxs[j[-1]]
                    self.txs[i].chans_to_pairs[dst] = Channel(dst, self.txs[i])
                    dst.chans_to_pairs[self.txs[i]] = self.txs[i].chans_to_pairs[dst]
                    self.txs[i].chans_to_pairs[dst].pow = j[1] * 1e3
                    self.txs[i].chans_to_pairs[dst].delay = j[2]
                    self.txs[i].chans_to_pairs[dst].delay_spread = j[3]
                    self.txs[i].chans_to_pairs[dst].dist = norm(self.txs[i].coords - dst.coords)
                    self.txs[i].chans_to_pairs[dst].chid = j[0]

                for j in self.txs[i].chans_to_pairs.keys():
                    self.dbcurs.execute(CHAN_PTH.format(self.txs[i].chans_to_pairs[j].chid))
                    d = self.dbcurs.fetchall()
                    d = sorted(d, key=lambda t: t[1], reverse=True)
                    for k in d[0:(self.npaths if d.__len__() >= self.npaths else d.__len__())]:
                        self.txs[i].chans_to_pairs[j].paths[k[0]] = Path()
                        self.txs[i].chans_to_pairs[j].paths[k[0]].pathid = k[0]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].chan = self.txs[i].chans_to_pairs[j]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].pow = k[1] * 1e3
                        self.txs[i].chans_to_pairs[j].paths[k[0]].fspl = k[7]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].phase = k[8]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].delay = k[2]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].aod = k[3]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].eod = k[4]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].aoa = k[5]
                        self.txs[i].chans_to_pairs[j].paths[k[0]].eoa = k[6]

        print('Success!')

        self.__unlink()

    def load_interactions(self, store: bool = True):
        print('Loading interactions...', end='', flush=True)

        self.__link()

        if hasattr(self, 'host') and self.threaded:
            txs_per_thread = self.txs.__len__() / self.nthreads
            rxs_per_thread = self.rxs.__len__() / self.nthreads

            thread_pool = cof.ThreadPoolExecutor(max_workers=self.nthreads)

            if txs_per_thread > rxs_per_thread:
                print('TXward...', end='', flush=True)
                txs = []
                txs_per_thread = int(txs_per_thread)
                for i in self.txs.keys():
                    if txs.__len__() % txs_per_thread == 0 and txs.__len__() != 0:
                        thread_pool.submit(_load_iters_txs, self, txs, self.host, self.user, self.pw,
                                               self.dbname, store)
                        txs = list()
                    txs.append(i)

                if txs.__len__() != 0 or txs_per_thread <= 1:
                    thread_pool.submit(_load_iters_txs, self, txs, self.host, self.user, self.pw, self.dbname, store)
            else:
                print('RXward...', end='', flush=True)
                rxs = []
                rxs_per_thread = int(rxs_per_thread)
                for i in self.rxs.keys():
                    if rxs.__len__() % rxs_per_thread == 0 and rxs.__len__() != 0:
                        thread_pool.submit(_load_iters_rxs, self, rxs, self.host, self.user, self.pw,
                                           self.dbname, store)
                        rxs = list()
                    rxs.append(i)

                if rxs.__len__() != 0 or rxs_per_thread <= 1:
                    thread_pool.submit(_load_iters_rxs, self, rxs, self.host, self.user, self.pw, self.dbname, store)

            thread_pool.shutdown()
        else:
            for i in self.txs.items():
                for j in i[1].chans_to_pairs.items():
                    self.dbcurs.execute(INTERS_SPEC_CHAN.format(j[1].chid))
                    inters = self.dbcurs.fetchall()
                    for k in j[1].paths.items():
                        dists = list()
                        precoords = i[1].coords
                        for l in inters:
                            if l[4] == k[1].pathid:
                                if store:
                                    intr = Interaction()
                                    intr.path = k[1]
                                    intr.coords = asarray([l[0], l[1], l[2]])
                                    intr.typ = l[3]
                                    dists.append(norm(precoords - intr.coords))
                                    precoords = intr.coords
                                else:
                                    intr = False
                                k[1].interactions.append(intr)
                        dists.append(norm(precoords - j[1].dest.coords))
                        k[1].length = sum(dists)
        print('Success!', flush=True)

        self.__unlink()

    def dump_paths(self,  txgrp: list = [-1], rxgrp: list = [-1], csvsav: bool = True, matsav: bool = True):
        for i in self.txs.items():
            if i[1].setid in txgrp or txgrp[0] == -1:
                for j in self.rxs.items():
                    aoa = list()
                    eoa = list()
                    aod = list()
                    eod = list()
                    phase = list()
                    power = list()
                    delay = list()
                    pathlen = list()
                    if j[1].setid in rxgrp or rxgrp[0] == -1:
                        if i[1].chan_to(j[1]):
                            for k in i[1].chans_to_pairs[j[1]].paths.items():
                                aoa.append(k[1].aoa)
                                eoa.append(k[1].eoa)
                                aod.append(k[1].aod)
                                eod.append(k[1].eod)
                                phase.append(k[1].phase)
                                power.append(k[1].pow)
                                delay.append(k[1].delay)
                                pathlen.append(k[1].length)
                        else:
                            pass

                    if matsav:
                        sio.savemat('Paths@[TX{}<->RX{}].mat'.format(i[0], j[0]), {'delay': delay, 'pow': l2db(power),
                                    'phase': phase, 'aoa': aoa, 'eoa': eoa, 'aod': aod, 'eod': eod, 'length': pathlen})

                    if csvsav:
                        file = open('Paths@[TX{}<->RX{}].csv'.format(i[0], j[0]), mode='w')
                        file.write('Delay [sec], Power [dBm], Phase, aoa [deg], eoa [deg], aod [deg], eod [deg], Length'
                                   ' [meters]\n')
                        for k in range(power.__len__()):
                            file.write('{},{},{},{},{},{},{},{}\n'.format(delay[k], l2db(power[k]), phase[k], aoa[k],
                                                                          eoa[k], aod[k], eod[k], pathlen[k]))
                        file.close()

    def load_all(self, dbname: str = ''):
        self.load_rxtx()
        self.load_paths()
        self.load_interactions()

    def end(self):
        self.dbcurs.close()
        self.dbconn.close()

    def __end__(self):
        self.end()



if __name__ == '__main__':
    DS = DataStorage(conf='dbconf.txt', dbname='Human_sitting_legsback_Sitting_sqlite')
    DS.load_rxtx()
    DS.load_paths()
    DS.load_interactions()
    from phys_path_procs import *
    check_data_nf(DS)
    exit()
