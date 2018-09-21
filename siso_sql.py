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


TX_EXTR = 'SELECT tx_id, x, y, z, tx_set_id FROM tx;'

RX_EXTR = 'SELECT rx_id, x, y, z, rx_set_id FROM rx;'

TX_PAIRS = 'SELECT * FROM (SELECT channel_utd.channel_id, channel_utd.received_power, ' \
           'channel_utd.mean_time_of_arrival, channel_utd.delay_spread ' \
           'FROM channel_utd WHERE channel_utd.channel_id IN ' \
           '(SELECT channel_id FROM channel WHERE tx_id = {})) utd ' \
           'JOIN ' \
           '(SELECT channel.channel_id ,channel.rx_id FROM channel WHERE tx_id = {}) chan ' \
           'ON utd.channel_id = chan.channel_id;'

TX_PAIRST = 'SELECT * FROM (SELECT channel_utd.channel_id, channel_utd.received_power, ' \
           'channel_utd.mean_time_of_arrival, channel_utd.delay_spread ' \
           'FROM channel_utd WHERE channel_utd.channel_id IN ' \
           '(SELECT channel_id FROM channel WHERE tx_id BETWEEN {} AND {})) utd ' \
           'JOIN ' \
           '(SELECT channel.channel_id ,channel.rx_id, channel.tx_id FROM channel WHERE tx_id BETWEEN {} AND {}) chan '\
           'ON utd.channel_id = chan.channel_id;'

RX_PAIRST = 'SELECT * FROM (SELECT channel_utd.channel_id, channel_utd.received_power, ' \
           'channel_utd.mean_time_of_arrival, channel_utd.delay_spread ' \
           'FROM channel_utd WHERE channel_utd.channel_id IN ' \
           '(SELECT channel_id FROM channel WHERE rx_id BETWEEN {} AND {})) utd ' \
           'JOIN ' \
           '(SELECT channel.channel_id ,channel.tx_id, channel.rx_id FROM channel WHERE rx_id BETWEEN {} AND {}) chan '\
           'ON utd.channel_id = chan.channel_id;'

CHAN_PTH = 'SELECT path_utd_id, received_power, time_of_arrival, departure_phi, departure_theta, arrival_phi,' \
           ' arrival_theta, freespace_path_loss, cir_phs FROM path_utd WHERE path_id IN (SELECT path_id FROM' \
           ' path WHERE channel_id = {});'

INTERS = 'SELECT * FROM interaction_type;'

INTERS_SPEC_CHAN = 'SELECT x, y, z, interaction_type_id, path_id FROM interaction WHERE path_id IN' \
                   '(SELECT path_id FROM path WHERE channel_id = {});'

INTERS_SPEC = 'SELECT x,y,z,interaction_type_id FROM interaction WHERE path_id = {};'
