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

import pairdata
import cir3d
import dirpaz
import asap_probs_extract

print('Loading data...')
data_storage = pairdata.data_stor()
data_storage.load_rxtx('/home/aleksei/Nextcloud/Documents/TTY/WORK/mmWave/Simulations/WI/Class@60GHz/TESTe/Class@60GHz.TESTe.sqlite')
data_storage.load_path()
print('Plotting 3D CIRs')
c3d = cir3d.cirs(data_storage)
c3d.draw(print=True)
print('Plotting reception patterns')
daz = dirpaz.rec_pat(data_storage)
daz.draw(print=True)
print('Printing distanced histogram')
asap = asap_probs_extract.distanced_hist_extractor(data_storage)
asap.build(typ='LOS')
asap.plot_hist()