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

import mysql.connector as msc
from mysql.connector.constants import ClientFlag
import sqlite3
import concurrent.futures as cof
from multiprocessing import cpu_count
from time import sleep
from sys import argv
from os import sep

max_waiting = 5 * cpu_count()
stepping = 1e3

WI_TABLES = {'path': [('path_id', 'INT PRIMARY KEY'), ('channel_id', 'INT'), ('foliage_distance', 'REAL')],
             'path_utd': [('path_id', 'INT'), ('path_utd_id', 'INT PRIMARY KEY'), ('tx_sub_antenna', 'INT'), ('rx_sub_antenna', 'INT'),
                          ('utd_instance_id', 'INT'), ('received_power', 'REAL'), ('time_of_arrival', 'REAL'),
                          ('departure_phi', 'REAL'), ('departure_theta', 'REAL'), ('arrival_phi', 'REAL'), ('arrival_theta', 'REAL'),
                          ('cmp_e_x_r', 'REAL'), ('cmp_e_y_r', 'REAL'), ('cmp_e_z_r', 'REAL'), ('cmp_e_x_i', 'REAL'), ('cmp_e_y_i', 'REAL'), ('cmp_e_z_i', 'REAL'),
                          ('cmp_h_x_r', 'REAL'), ('cmp_h_y_r', 'REAL'), ('cmp_h_z_r', 'REAL'), ('cmp_h_x_i', 'REAL'), ('cmp_h_y_i', 'REAL'), ('cmp_h_z_i', 'REAL'),
                          ('freespace_path_loss', 'REAL'), ('freespace_path_loss_woa', 'REAL'), ('e_theta_r', 'REAL'), ('e_theta_i', 'REAL'),
                          ('e_phi_r', 'REAL'), ('e_phi_i', 'REAL'), ('cir_phs', 'REAL'), ('cmp_volt_r', 'REAL'), ('cmp_volt_i', 'REAL')],
             'channel': [('channel_id', 'INT PRIMARY KEY'), ('tx_id', 'INT'), ('rx_id', 'INT')],
             'channel_utd': [('channel_utd_id', 'INT PRIMARY KEY'), ('channel_id', 'INT'), ('utd_instance_id', 'INT'), ('received_power', 'REAL'),
                             ('mean_time_of_arrival', 'REAL'), ('delay_spread', 'REAL')],
             'rx': [('rx_id', 'INT PRIMARY KEY'), ('rx_set_id', 'INT'), ('x', 'REAL'), ('y', 'REAL'), ('z', 'REAL')],
             'tx': [('tx_id', 'INT PRIMARY KEY'), ('tx_set_id', 'INT'), ('x', 'REAL'), ('y', 'REAL'), ('z', 'REAL')],
             'diffraction_edge': [('edge_id', 'INT PRIMARY KEY'), ('vertex_a', 'INT'), ('vertex_b', 'INT'), ('edgefacet_a', 'INT'), ('edgefacet_b', 'INT'), ('angle', 'REAL')],
             'interaction': [('interaction_id', 'INT PRIMARY KEY'), ('path_id', 'INT'), ('interaction_type_id', 'INT'), ('object_id', 'INT'), ('sub_id', 'INT'), ('x', 'REAL'), ('y', 'REAL'), ('z', 'REAL'),
                             ('has_location', 'INT')],
             'interaction_type': [('interaction_type_id', 'INT PRIMARY KEY'), ('description', 'CHAR(64)')],
             'rx_set': [('rx_set_id', 'INT PRIMARY KEY'), ('txrx_set_type_id', 'INT'), ('spacing1', 'REAL'), ('spacing2', 'REAL'), ('spacing3', 'REAL')],
             'rx_metadata': [('rx_metadata_id', 'INT PRIMARY KEY'), ('rx_id', 'INT'), ('utd_instance_id', 'INT'), ('max_gain', 'REAL')],
             'tx_metadata': [('tx_metadata_id', 'INT'), ('tx_id', 'INT'), ('utd_instance_id', 'INT'), ('max_gain', 'REAL'), ('power', 'REAL')],
             'txrx_set_type': [('txrx_set_type_id', 'INT PRIMARY KEY'), ('description', 'TEXT')],
             'utd_instance': [('utd_instance_id', 'INT PRIMARY KEY')],
             'utd_instance_param': [('utd_instance_id', 'INT NOT NULL'), ('utd_instance_param_id', 'INT NOT NULL'), ('element_id', 'INT NOT NULL'), ('parameter', 'REAL')],
             'scene_origin': [('latitude', 'REAL'), ('longitude', 'REAL'), ('altitude', 'REAL')]}

conf = open('dbconf.txt')

host = conf.readline()
user = conf.readline()
pw = conf.readline()

conf.close()

sqlconn = msc.connect(host=host.strip('\n'), user=user.strip('\n'), password=pw.strip('\n'), client_flags=[ClientFlag.SSL])
sqlcurs = sqlconn.cursor()


if argv.__len__() == 1:
    dbf = input('Type in DB path: ')
    dbn = input('Type in database name: ')
else:
    dbf = ' '.join(argv[1:argv.__len__()]).strip('"').strip("'")
    dbn = dbf.split(sep)[-1]

dbn = dbn.replace('.', '_').replace('@','at').replace(' ', '_')

sqlcurs.execute('CREATE DATABASE {};'.format(dbn))
sqlcurs.execute('USE {};'.format(dbn))

# Create table on MySQL server
for i in WI_TABLES.items():
    # Compile request
    if i[0] in ['path', 'path_utd', 'channel', 'channel_utd', 'rx', 'tx', 'diffraction_edge', 'interaction', 'rx_set',
                'rx_metadata', 'tx_metadata', 'txrx_set_type', 'utd_instance', 'utd_instance_param', 'scene_origin']:
        rstr = 'CREATE TABLE {}('.format(i[0])
        for j in i[1]:
            rstr = rstr + j[0] + ' ' + j[1] + (',\n' if i[1][-1] != j else '\n')

        rstr = rstr + ');'
        sqlcurs.execute(rstr)

print(dbf)
sleep(5)

sqlicon = sqlite3.connect(dbf)
sqlicon.row_factory = sqlite3.Row
sqlicurs = sqlicon.cursor()

waiting = 0


def sqlins(dbn, data, host, user, pw):
    global waiting
    waiting += 1
    conn = msc.connect(host=host.strip('\n'), user=user.strip('\n'), password=pw.strip('\n'))
    curs = conn.cursor()
    curs.execute('USE {};'.format(dbn))
    curs.execute(data)
    conn.commit()
    conn.close()
    waiting -= 1


TPE = cof.ThreadPoolExecutor(max_workers=max_waiting)

for i in WI_TABLES.items():
    if i[0] in ['path', 'path_utd', 'channel', 'channel_utd', 'rx', 'tx', 'diffraction_edge', 'interaction', 'rx_set',
                'rx_metadata', 'tx_metadata', 'txrx_set_type', 'utd_instance', 'utd_instance_param', 'scene_origin']:
        # Construct request
        req = "SELECT "
        for j in i[1]:
            req = req + j[0] + (', ' if i[1][-1] != j else ' FROM {}'.format(i[0]))

        comm = 0

        key_prep = False

        print('Writing {}'.format(i[0]))

        myreq = []

        for j in sqlicurs.execute(req):
            if comm % stepping == 0:
                myreq = ["INSERT INTO {} ".format(i[0])]
                myreq.append('(')
                myreq.append(','.join(j.keys()))
                myreq.append(') VALUES ')

            myreq.append('(')

            te = 1
            for k in j.keys():
                myreq.append(str(j[k]))
                myreq.append((',' if te < j.keys().__len__() else ')'))
                te += 1

            if (comm + 1) % stepping == 0:
                myreq.append(';')
                # Avoid excessive memory consumption
                while waiting >= max_waiting:
                    sleep(0.1)

                TPE.submit(sqlins, dbn, ''.join(myreq), host, user, pw)
                print('Writing to {} at row {}'.format(i[0], comm))
            else:
                myreq.append(',')

            comm += 1

        if (comm + 1) % stepping != 0:
            myreq[-1] = ';'
            while waiting >= max_waiting:
                sleep(0.1)
            TPE.submit(sqlins, dbn, ''.join(myreq), host, user, pw)
            print('Writing to {} at row {}'.format(i[0], comm))

        print('Finished writing {}'.format(i[0]))

TPE.shutdown(wait=True)

sqlcurs.execute('CREATE INDEX channel_tx_rx_index ON channel(tx_id,rx_id);')
sqlcurs.execute('CREATE INDEX channel_utd_channel_index ON channel_utd(channel_id);')
sqlcurs.execute('CREATE INDEX interaction_path_index ON interaction(path_id);')
sqlcurs.execute('CREATE INDEX path_channel_index ON path(channel_id);')
sqlcurs.execute('CREATE INDEX path_utd_path_index ON path_utd(path_id);')

sqlconn.commit()
sqlconn.close()
print('ALL DONE, CONNECTIONS TERMINATED')
input('')
exit()
