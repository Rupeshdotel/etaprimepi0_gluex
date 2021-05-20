#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 14:39:09 2021

@author: rupeshdotel
"""

import ROOT as R
from root_numpy import tree2array
import LT.box as B
#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree__030274_cons.root")
intree = rfile.Get('qfactortree')
dc = tree2array(intree)


#%%
mpi0c = dc['mpi0']
mpi0mc = dc['mpi0m']

metac = dc['meta']
metamc = dc['metam']


h_pi0c = B.histo(mpi0c, bins = 30, title = "Inavriant mass gg kinfit variables")
h_pi0mc = B.histo(mpi0mc, bins = 30, range = [0.07, 0.22], title = "Inavriant mass gg Measured variables")

h_etac = B.histo(metac, bins = 30, title = "Inavriant mass gg kinfit variables")
h_etamc = B.histo(metamc, bins = 30, range = [0.30, 0.80], title = "Inavriant mass gg Measured variables")



#%%
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree__030274_uncons.root")
intree = rfile.Get('qfactortree')
d = tree2array(intree)

#%%
mpi0 = d['mpi0']
mpi0m = d['mpi0m']

meta = d['meta']
metam = d['metam']


h_pi0 = B.histo(mpi0, bins = 30, range = [0.00, 0.40], title = "Inavriant mass gg kinfit variables")
h_pi0m = B.histo(mpi0m, bins = 30, range = [0.07, 0.22], title = "Inavriant mass gg Measured variables")

h_eta = B.histo(meta, bins = 30, title = "Inavriant mass gg kinfit variables")
h_etam = B.histo(metam, bins = 30, range = [0.30, 0.80], title = "Inavriant mass gg Measured variables")








