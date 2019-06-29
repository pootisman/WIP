import sqlite3
import mysql.connector as msqlc
from mysql.connector.constants import ClientFlag
from siso_sql import *
from numpy import zeros, asarray
from numpy.linalg import norm


class Node:
    def __init__(self, typ: str):
        self.chans_to_pairs = None
        self.node_id = 0
        self.coords = zeros([3])
        self.rot = zeros([3])
        self.setid = 0
        self.type = typ

        if typ == 'TX':
            self.txpow = 0.0
        else:
            self.rxpow = 0.0

    def chan_to(self, dest):
        if self.type == 'TX':
            for i in self.chans_to_pairs.keys():
                if self.chans_to_pairs[i].dest == dest:
                    return self.chans_to_pairs[i]
        else:
            for i in self.chans_to_pairs.keys():
                if self.chans_to_pairs[i].src == dest:
                    return self.chans_to_pairs[i]
        return None

    def __repr__(self):
        return '{} node {}'.format(self.type, self.node_id)


class Channel:
    def __init__(self, dest: Node = None, src: Node = None):
        self.paths = dict()
        self.dest = dest
        self.src = src
        self.pow = 0.0
        self.delay = 0.0
        self.delay_spread = 0.0
        self.dist = 0.0
        self.chid = 0
        self.clusters = dict()

    def __repr__(self):
        return 'Radio channel'

    def __str__(self):
        return 'Sumpow = {}, delay = {}, delay_spread = {}, Ncl = {}, {} -> {}'.\
            format(self.pow, self.delay, self.delay_spread, self.clusters.__len__(), self.src.node_id, self.dest.node_id)


class Path:
    def __init__(self):
        self.pathid = 0
        self.pow = 0.0
        self.phase = 0.0
        self.delay = 0.0
        self.len = 0.0
        self.interactions = list()
        self.aoa = 0.0
        self.eoa = 0.0
        self.aod = 0.0
        self.eod = 0.0
        self.fspl = 0.0
        self.length = 0.0
        self.chan = None
        self.cluster = None
        self.near_field_failed = False

    def __repr__(self):
        return 'Propagation path'

    def __str__(self):
        return 'Inters = {}, pow = {} mW, len = {} m, {} [Az:{},El:{}]-->[Az:{},El:{}] {}'.\
            format(self.interactions.__len__(), self.pow, self.len, self.chan.src.node_id, self.aod, self.eod,
                   self.aoa, self.eoa, self.chan.dest.node_id)


class Interaction:
    def __init__(self):
        self.typ = 'TX'
        self.coords = zeros([3])
        self.path = None

    def __repr__(self):
        print('Interaction')

    def __str__(self):
        print('{}@{}'.format(self.typ, self.coords))


class RXConnector:
    def __init__(self, dbcursor: None, master = None):
        self.dbcurs = dbcursor
        self.master = master

    def items(self, grprange: list = [-1], rxrange: list = [-1]):
        output = []

        reqstr = list(RX_EXTR)

        if grprange != [-1] and rxrange == [-1]:
            reqstr[-1] = ''
            reqstr = ''.join(reqstr) + ' WHERE rx_set_id IN {};'.format('({})'.format(','.join('%d' % i for i in grprange)))
        elif grprange == [-1] and rxrange != [-1]:
            reqstr[-1] = ''
            reqstr = ''.join(reqstr) + ' WHERE rx_id IN {};'.format('({})'.format(','.join('%d' % i for i in rxrange)))
        elif grprange != [-1] and rxrange != [-1]:
            reqstr[-1] = ''
            reqstr = ''.join(reqstr) + ' WHERE rx_id IN {} AND rx_set_id IN {};'.format('({})'.format(','.join('%d' % i for i in rxrange)),
                                                                               '({})'.format(','.join('%d' % i for i in grprange)))
        else:
            reqstr = ''.join(reqstr)

        self.dbcurs.execute(reqstr, multi=True)
        data = self.dbcurs.fetchall()

        for i in data:
            n = Node(typ='RX')
            n.chans_to_pairs = ChannelConnector(dbcursor=self.dbcurs)
            n.chans_to_pairs.origin = n
            n.node_id = i[0]
            n.coords = asarray([i[1], i[2], i[3]])
            n.setid = i[4]
            output.append((i[0], n))

        print(output)

        return output

    def __getitem__(self, key):
        reqstr = list(RX_EXTR)
        reqstr[-1] = ''
        reqstr = ''.join(reqstr) + ' WHERE rx_id = {};'.format(key)

        self.dbcurs.execute(reqstr, multi=True)
        data = self.dbcurs.fetchall()
        if len(data) == 0:
            raise KeyError
        n = Node(typ='RX')
        data = data[0]
        n.node_id = data[0]
        n.coords = asarray([data[1], data[2], data[3]])
        n.chans_to_pairs = ChannelConnector(dbcursor=self.dbcurs, origin=n, master=self)
        n.setid = data[4]

        print(n)

        return n


class TXConnector:
    def __init__(self, dbcursor: None, master = None):
        self.master = master
        self.dbcurs = dbcursor

    def items(self, grprange: list = [-1], txrange: list = [-1]):
        output = []

        reqstr = list(TX_EXTR)

        if grprange != [-1] and txrange == [-1]:
            reqstr[-1] = ''
            reqstr = ''.join(reqstr) + ' WHERE rx_set_id IN {};'.format('({})'.format(','.join('%d' % i for i in grprange)))
        elif grprange == [-1] and txrange != [-1]:
            reqstr[-1] = ''
            reqstr = ''.join(reqstr) + ' WHERE rx_id IN {};'.format('({})'.format(','.join('%d' % i for i in txrange)))
        elif grprange != [-1] and txrange != [-1]:
            reqstr[-1] = ''
            reqstr = ''.join(reqstr) + ' WHERE rx_id IN {} AND rx_set_id IN {};'.format('({})'.format(','.join('%d' % i for i in txrange)),
                                                                               '({})'.format(','.join('%d' % i for i in grprange)))
        else:
            reqstr = ''.join(reqstr)

        self.dbcurs.execute(reqstr, multi=True)
        data = self.dbcurs.fetchall()

        for i in data:
            n = Node(typ='TX')
            n.chans_to_pairs = ChannelConnector(dbcursor=self.dbcurs, origin=n, master=self)
            n.chans_to_pairs.origin = n
            n.node_id = i[0]
            n.coords = asarray([i[1], i[2], i[3]])
            n.setid = i[4]
            output.append((i[0], n))

        print(output)

        return output

    def __getitem__(self, key):
        reqstr = list(TX_EXTR)
        reqstr[-1] = ''
        reqstr = ''.join(reqstr) + ' WHERE tx_id = {};'.format(key)

        self.dbcurs.execute(reqstr, multi=True)
        data = self.dbcurs.fetchall()
        if len(data) == 0:
            raise KeyError
        n = Node(typ='TX')
        data = data[0]
        n.node_id = data[0]
        n.coords = asarray([data[1], data[2], data[3]])
        n.chans_to_pairs = ChannelConnector(dbcursor=self.dbcurs, origin=n, master=self)
        n.setid = data[4]

        print(n)

        return n


class ChannelConnector:
    def __init__(self, dbcursor: None, origin: Node = None, destination: Node = None, master = None):
        self.master = master
        self.origin = origin
        self.dest = destination
        self.dbcurs = dbcursor

    def items(self, destrange: list = [-1], destgrprange: list = [-1]):
        output = []
        additional_filter = ''

        if self.origin.type == 'TX':
            reqstr = TX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND rx_id IN ({})'.format(','.join('%d' % i for i in destrange))
            if destgrprange != [-1]:
                additional_filter = additional_filter + ' AND rx_id IN (SELECT rx_id FROM rx WHERE rx_set_id IN ({}))'.format(','.join('%d' % i for i in destgrprange))
        else:
            reqstr = RX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND tx_id IN ({})'.format('%d' % i for i in destrange)

            if destgrprange != [-1]:
                additional_filter = additional_filter + ' AND tx_id IN (SELECT tx_id FROM tx WHERE tx_set_id IN ({}))'.format(','.join('%d' % i for i in destgrprange))

        reqstr = reqstr.format(self.origin.node_id, additional_filter, self.origin.node_id)

        self.dbcurs.execute(reqstr, multi=True)

        data = self.dbcurs.fetchall()

        for i in data:
            chan = Channel()
            chan.delay_spread = i[3]
            chan.delay = i[2]
            chan.pow = i[1]
            chan.chid = i[0]
            chan.paths = PathConnector(dbcursor=self.dbcurs, master=self, partof=chan)
            if self.origin.type == 'TX':
                chan.src = self.origin
                chan.dest = self.master.master.rxs[i[6]]
                chan.dist = norm(self.origin.coords - self.master.master.rxs[i[6]].coords)
            else:
                chan.dest = self.origin
                chan.src = self.master.master.txs[i[5]]
                chan.dist = norm(self.origin.coords - self.master.master.txs[i[5]].coords)
            output.append((chan.dest, chan))

        print(output)

        return output

    def keys(self, destrange: list = [-1]):
        output = []
        additional_filter = ''

        if self.origin.type == 'TX':
            reqstr = TX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND rx_id IN ({}) '.format(','.join('%d' % i for i in destrange))
        else:
            reqstr = RX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND tx_id IN ({}) '.format(','.join('%d' % i for i in destrange))

        reqstr = reqstr.format(self.origin.node_id, additional_filter, self.origin.node_id)

        self.dbcurs.execute(reqstr, multi=True)

        data = self.dbcurs.fetchall()

        for i in data:
            output.append(i[0])

        print(output)

        return output

    def __getitem__(self, key: int = 0):
        chan = Channel()

        if self.origin.type == 'TX':
            additional_filter = ' AND rx_id = {}'.format(key)
            reqstr = TX_CHAN_CONNECTOR
        else:
            additional_filter = ' AND tx_id = {}'.format(key)
            reqstr = RX_CHAN_CONNECTOR

        reqstr = reqstr.format(self.origin.node_id, additional_filter, self.origin.node_id)

        self.dbcurs.execute(reqstr, multi=True)

        data = self.dbcurs.fetchall()

        if len(data) == 0:
            raise KeyError

        data = data[0]

        chan.delay_spread = data[3]
        chan.delay = data[2]
        chan.pow = data[1]
        chan.chid = data[0]
        chan.paths = PathConnector(dbcursor=self.dbcurs, master=self, partof=chan)
        if self.origin.type == 'TX':
            chan.src = self.origin
            chan.dest = self.master.master.rxs[key]
            chan.dist = norm(self.origin.coords - self.master.master.rxs[key].coords)
        else:
            chan.dest = self.origin
            chan.src = self.master.master.txs[key]
            chan.dist = norm(self.origin.coords - self.master.master.txs[key].coords)

        print(chan)

        return chan


class PathConnector:
    def __init__(self, dbcursor: None, partof: Channel = None, master = None):
        self.master = master
        self.partof = partof
        self.dbcurs = dbcursor

    def items(self):
        output = []

        reqstr = CHAN_PTH.format(self.partof.chid)

        self.dbcurs.execute(reqstr, multi=True)
        data = self.dbcurs.fetchall()

        for i in data:
            p = Path()
            p.pathid = i[0]
            p.pow = i[1]
            p.delay = i[2]
            p.aod = i[3]
            p.eod = i[4]
            p.aoa = i[5]
            p.eoa = i[6]
            p.fspl = i[7]
            p.phase = i[8]
            output.append((p.pathid, p))

        print(output)

        return output

    def __getitem__(self, key):
        reqstr = CHAN_PTH.format(self.partof.chid)

        self.dbcurs.execute(reqstr)
        data = self.dbcurs.fetchall()

        data = data[0]

        p = Path()
        p.pathid = data[0]
        p.pow = data[1]
        p.delay = data[2]
        p.aod = data[3]
        p.eod = data[4]
        p.aoa = data[5]
        p.eoa = data[6]
        p.fspl = data[7]
        p.phase = data[8]

        print(p)

        return p