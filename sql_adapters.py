import sqlite3
import mysql.connector as msqlc
from pairdata import *
from mysql.connector.constants import ClientFlag
from siso_sql import *


class RXConnector:
    def __init__(self, db: str = 'test', host: str = '', user: str = 'user', pw: str = 'passwd'):
        self.db = db
        self.host = host
        self.user = user
        self.user = user
        self.pw = pw

        if host == '':
            self.dbconn = sqlite3.connect(db)
            self.dbcurs = self.dbconn.cursor()
        else:
            self.dbconn = msqlc.connect(host=host, user=user, password=pw,
                                        client_flags=[ClientFlag.SSL], database=db)
            self.dbcurs = self.dbconn.cursor()
            self.dbcurs.execute('USE {};'.format(db))

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

        print(reqstr)
        self.dbcurs.execute(reqstr)
        data = self.dbcurs.fetchall()

        for i in data:
            n = Node(typ='RX')
            n.chans_to_pairs = ChannelConnector(db = self.db, host = self.host,
                                                user = self.user, pw = self.pw,
                                                node=n)
            n.node_id = i[0]
            n.x = i[1]
            n.y = i[2]
            n.z = i[3]
            n.setid = i[4]
            output.append((i[0],n))

        return output


class TXConnector:
    def __init__(self, db: str = 'test', host: str = '', user: str = 'user', pw: str = 'passwd'):
        self.db = db
        self.host = host
        self.user = user
        self.user = user
        self.pw = pw

        if host == '':
            self.dbconn = sqlite3.connect(db)
            self.dbcurs = self.dbconn.cursor()
        else:
            self.dbconn = msqlc.connect(host=host, user=user, password=pw,
                                        client_flags=[ClientFlag.SSL], database=db)
            self.dbcurs = self.dbconn.cursor()
            self.dbcurs.execute('USE {};'.format(db))

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

        print(reqstr)
        self.dbcurs.execute(reqstr)
        data = self.dbcurs.fetchall()

        for i in data:
            n = Node(typ='TX')
            n.chans_to_pairs = ChannelConnector(db = self.db, host = self.host,
                                                user = self.user, pw = self.pw,
                                                node=n)
            n.node_id = i[0]
            n.x = i[1]
            n.y = i[2]
            n.z = i[3]
            n.setid = i[4]
            output.append((i[0],n))

        return output


class ChannelConnector:
    def __init__(self, db: str = 'test', host: str = '', user: str = 'user', pw: str = 'passwd',
                 node: Node = None):
        self.db = db
        self.host = host
        self.user = user
        self.user = user
        self.pw = pw
        self.origin = node

        if host == '':
            self.dbconn = sqlite3.connect(db)
            self.dbcurs = self.dbconn.cursor()
        else:
            self.dbconn = msqlc.connect(host=host, user=user, password=pw,
                                        client_flags=[ClientFlag.SSL], database=db)
            self.dbcurs = self.dbconn.cursor()
            self.dbcurs.execute('USE {};'.format(db))

    def items(self, destrange: list = [-1]):
        output = []
        additional_filter = ''

        if self.origin.type == 'TX':
            reqstr = TX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND rx_id IN ({})'.format('%d' % i for i in destrange)
        else:
            reqstr = RX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND tx_id IN ({})'.format('%d' % i for i in destrange)

        reqstr.format(self.origin.node_id, additional_filter, self.origin.node_id)

        self.dbcurs.execute(reqstr)

        data = self.dbcurs.fetchall()

        for i in data:
            c = Channel()
            c.paths =PathConnector(db = self.db, host = self.host,
                                   user = self.user, pw = self.pw,
                                   channel=c)
            output.append(i[0], c)

        return output

    def keys(self, destrange: list = [-1]):
        output = []
        additional_filter = ''

        if self.origin.type == 'TX':
            reqstr = TX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND rx_id IN ({}) '.format('%d' % i for i in destrange)
        else:
            reqstr = RX_CHAN_CONNECTOR

            if destrange != [-1]:
                additional_filter = ' AND tx_id IN ({}) '.format('%d' % i for i in destrange)

        reqstr.format(self.origin.node_id, additional_filter, self.origin.node_id)

        self.dbcurs.execute(reqstr)

        data = self.dbcurs.fetchall()

        for i in data:
            output.append(i[0])

        return output

    def __getitem__(self, dest: int = 0):
        chan = Channel()

        additional_filter = ' AND tx_id IN ({}) '.format(dest)

        if self.origin.type == 'TX':
            reqstr = TX_CHAN_CONNECTOR
        else:
            reqstr = RX_CHAN_CONNECTOR

        reqstr.format(self.origin.node_id, additional_filter, self.origin.node_id)

        self.dbcurs.execute(reqstr)

        data = self.dbcurs.fetchall()

        for i in data:
            chan.delay = i[2]
            chan.pow = i[1]
            chan.chid = i[0]
            chan.paths = PathConnector(db = self.db, host = self.host,
                                    user = self.user, pw = self.pw,
                                    channel=chan)

        return chan


class PathConnector:
    def __init__(self, db: str = 'test', host: str = '', user: str = 'user', pw: str = 'passwd', channel: Channel = None):
        self.db = db
        self.host = host
        self.user = user
        self.user = user
        self.pw = pw

        if host == '':
            self.dbconn = sqlite3.connect(db)
            self.dbcurs = self.dbconn.cursor()
        else:
            self.dbconn = msqlc.connect(host=host, user=user, password=pw,
                                        client_flags=[ClientFlag.SSL], database=db)
            self.dbcurs = self.dbconn.cursor()
            self.dbcurs.execute('USE {};'.format(db))

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

        print(reqstr)
        self.dbcurs.execute(reqstr)
        data = self.dbcurs.fetchall()

        for i in data:
            n = Node(typ='TX')
            n.chans_to_pairs = ChannelConnector(db = self.db, host = self.host,
                                                user = self.user, pw = self.pw,
                                                node=n)
            n.node_id = i[0]
            n.x = i[1]
            n.y = i[2]
            n.z = i[3]
            n.setid = i[4]
            output.append((i[0],n))

        return output