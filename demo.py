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

from pairdata import DataStorage
from cir import CIR
from circollator import CIRCollator
from plpow import PLPlot
from auxfun import enable_latex
from asap_probs_extract import DistancedHistExtractor
from delay_spread import DelaySpreadPlot

gencir = True
gencircoll = True
genreg = False
disthis = False
ds = False

print('Loading data fleece...')
DS = DataStorage('dbconf.txt')
DS.load_rxtx('Human_sitting_legsback_Standing_fleece_sqlite')
DS.load_paths(npaths=250)
DS.load_interactions(store=True)
print('Loading data cotton...')
#DC = DataStorage('dbconf.txt')
#DC.load_rxtx('Human_sitting_legsback_Sitting_cotton_sqlite')
#DC.load_paths(npaths=250)
#DC.load_interactions(store=True)
print('Loading data Leather...')
DL = DataStorage('dbconf.txt')
DL.load_rxtx('Human_sitting_legsback_Standing_leather_sqlite')
DL.load_paths(npaths=250)
DL.load_interactions(store=True)
print('Loading data naked...')
DN = DataStorage('dbconf.txt')
DN.load_rxtx('Standing_naked_sqlite')
DN.load_paths(npaths=250)
DN.load_interactions(store=True)


rxgrp = 2

enable_latex()

if gencir:
    print('Plotting 3D CIRs')
    c3ds = CIR(DS)
    c3ds.export(cmap='Blues', rxgrp=rxgrp, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=1, title='Fleece ')
    #c3dc = CIR(DC)
    #c3dc.export(rxgrp=rxgrp, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=2, title='Cotton ')
    c3dl = CIR(DL)
    c3dl.export(rxgrp=rxgrp, mkpng=False, show=False, zmin=-110, zmax=-40, fidbase=3, title='Leather ')
    c3dn = CIR(DN)
    c3dn.export(cmap='Blues', rxgrp=rxgrp, mkpng=False, show=True, zmin=-110, zmax=-40, fidbase=4, title='Naked ')

    if gencircoll:
        print('Collating CIRs')
        ccl = CIRCollator()

        ccl + [c3dn, c3ds, c3dl]

        ccl.export_collated(show=True, idxs=[0, 1], csq=True, xlabel=False, ylabel=True, title_draw=False,
                            figsize=(5, 6), csqloc=2)
        ccl.export_collated(show=True, idxs=[0, 2], csq=True, xlabel=False, ylabel=False, title_draw=False,
                            figsize=(5, 6))
    #ccl.export_collated(show=True, idxs=[0, 2], csq=True, xlabel=False, ylabel=False, title_draw=False, figsize=(5, 6),
    #                    csqloc=2)
    #ccl.export_collated(show=True, idxs=[0, 2], csq=True, xlabel=False, ylabel=False, title_draw=False, figsize=(5, 6),
    #                    csqloc=2)


if genreg:
    print('Trying to build regression')
    plp = PLPlot(source=DN)
    print(plp.regr_comp(typ='NLOS-1', rxgrp=[2, 4, 5], threshold=-250, nff=False))
    plp.export(csvsav=True, matsav=True)

if disthis:
    print('Printing distanced histogram')
    asap = DistancedHistExtractor(DN, nffilt=False, thrs=-150, range=(0.1, 1.0), minbins=0.05)
    asap.build_trans(typ='LOS->LOS')
    asap.plot_hist()

if ds:
    print('Printing delaye spread')
    dspn = DelaySpreadPlot(DN)
    dspn.plot_groups(rxgrp=[2, 4, 5], overlay=True, title='Naked ', ymin=0.0, ymax=1.0, matsav=True,
                     rx_name_map={2: "Belt", 4: "Up", 5: "Diag"})
    dsps = DelaySpreadPlot(DS)
    dsps.plot_groups(rxgrp=[2, 4, 5], overlay=True, title='Fleece ', ymin=0.0, ymax=1.0, matsav=True,
                     rx_name_map={2: "Belt", 4: "Up", 5: "Diag"})
    dspl = DelaySpreadPlot(DL)
    dspl.plot_groups(rxgrp=[2, 4, 5], overlay=True, title='Leather ', ymin=0.0, ymax=1.0, matsav=True,
                     rx_name_map={2: "Belt", 4: "Up", 5: "Diag"})
    dspc = DelaySpreadPlot(DC)
    dspc.plot_groups(rxgrp=[2, 4, 5], overlay=True, title='Cotton ', ymin=0.0, ymax=1.0, matsav=True,
                     rx_name_map={2: "Belt", 4: "Up", 5: "Diag"})
