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
import cir
import circollator
import asap_probs_extract
import matplotlib.pyplot as mpl

print('Loading data fleece...')
DS = pairdata.data_stor('dbconf.txt')
DS.load_rxtx('Human_sitting_legsback_Sitting_fleece_sqlite')
DS.load_paths(npaths=250)
DS.load_interactions(store=True)
#print('Loading data cotton...')
#DC = pairdata.data_stor('dbconf.txt')
#DC.load_rxtx('Human_sitting_legsback_Sitting_cotton_sqlite')
#DC.load_paths(npaths=250)
#DC.load_interactions(store=True)
print('Loading data Leather...')
DL = pairdata.data_stor('dbconf.txt')
DL.load_rxtx('Human_sitting_legsback_Sitting_Leather_sqlite')
DL.load_paths(npaths=250)
DL.load_interactions(store=True)
print('Loading data naked...')
DN = pairdata.data_stor('dbconf.txt')
DN.load_rxtx('Human_sitting_legsback_Sitting_sqlite')
DN.load_paths(npaths=250)
DN.load_interactions(store=True)

print('Plotting 3D CIRs')
c3ds = cir.cirs(DS)
c3ds.export(rxgrp=6, mkpng=True, show=False, zmin=-140, zmax=-40, fidbase=1, title='Fleece ')
#c3dc = cir.cirs(DC)
#c3dc.export(rxgrp=6, mkpng=True, show=False, zmin=-140, zmax=-40, fidbase=2, title='Cotton ')
c3dl = cir.cirs(DL)
c3dl.export(rxgrp=6, mkpng=True, show=False, zmin=-140, zmax=-40, fidbase=2, title='Leather ')
c3dn = cir.cirs(DN)
c3dn.export(rxgrp=6, mkpng=True, show=False, zmin=-140, zmax=-40, fidbase=3, title='Naked ')

ccl = circollator.circollator()

ccl + [c3dn, c3ds]

ccl.export_collated(show=True)

#print('Printing distanced histogram')
#asap = asap_probs_extract.distanced_hist_extractor(DS)
#asap.build(typ='LOS')
#asap.plot_hist()