#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 11:20:44 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import ROOT as R
from root_numpy import tree2array
import matplotlib.pyplot as plt


#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/Tree_etapr_gpippim_0303.root")
intree = rfile.Get('Tree_etapr_gpippim')
d = tree2array(intree)

#%%

mm2m = d['mm2m']
mpippim = d['mpippim']
mgpippim = d['mgpippim']
mpipp = d['mpipp']
mpimp = d['mpimp']
mpippimp = d['mpippimp']

mant = d['mant']

#%%

plt.figure();
h = B.histo(mpippimp, bins = 50)
h.plot_exp()
plt.xlabel("")
plt.title("")
