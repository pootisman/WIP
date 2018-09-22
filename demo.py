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

import pairdata
import cir
import circollator
from auxfun import enable_latex

print('Loading data fleece...')
DS = pairdata.DataStorage('dbconf.txt')
DS.load_rxtx('Human_sitting_legsback_Sitting_fleece_sqlite')
DS.load_paths(npaths=250)
DS.load_interactions(store=True)
#print('Loading data cotton...')
#DC = pairdata.data_stor('dbconf.txt')
#DC.load_rxtx('Human_sitting_legsback_Sitting_cotton_sqlite')
#DC.load_paths(npaths=250)
#DC.load_interactions(store=True)
print('Loading data Leather...')
DL = pairdata.DataStorage('dbconf.txt')
DL.load_rxtx('Human_sitting_legsback_Sitting_Leather_sqlite')
DL.load_paths(npaths=250)
DL.load_interactions(store=True)
print('Loading data naked...')
DN = pairdata.DataStorage('dbconf.txt')
DN.load_rxtx('Human_sitting_legsback_Sitting_sqlite')
DN.load_paths(npaths=250)
DN.load_interactions(store=True)


rxgrp = 5

enable_latex()
print('Plotting 3D CIRs')
c3ds = cir.CIR(DS)
c3ds.export(rxgrp=rxgrp, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=1, title='Fleece ')
#c3dc = cir.cirs(DC)
#c3dc.export(rxgrp=2, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=2, title='Cotton ')
c3dl = cir.CIR(DL)
c3dl.export(rxgrp=rxgrp, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=2, title='Leather ')
c3dn = cir.CIR(DN)
c3dn.export(cmap='Blues', rxgrp=rxgrp, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=3,
            title='Naked ')

ccl = circollator.CIRCollator()

ccl + [c3dn, c3ds, c3dl]

ccl.export_collated(show=True, idxs=[0, 1], csq=True, xlabel=False, ylabel=True, title_draw=False, figsize=(5, 6),
                    csqloc=2)
ccl.export_collated(show=True, idxs=[0, 2], csq=True, xlabel=False, ylabel=False, title_draw=False, figsize=(5, 6),
                    csqloc=2)
#ccl.export_collated(show=True, idxs=[0, 2])
#ccl.export_collated(show=True, idxs=[0, 3])


#print('Printing distanced histogram')
#asap = asap_probs_extract.distanced_hist_extractor(DS)
#asap.build(typ='LOS')
#asap.plot_hist()