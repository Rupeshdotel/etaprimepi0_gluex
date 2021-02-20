#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 20:13:21 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import ROOT as R
from root_numpy import tree2array
import matplotlib.pyplot as plt

#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/gen_mctree.root")
intree = rfile.Get('gen_mctree')

d = tree2array(intree)

   
#%%   
beam_energy = d['beam_energy']
mproton = d['mproton']
metap = d['metap']
mpi0 = d['mpi0']
metappi0 = d['metaprimepi0']
mpi0p = d['mpi0p']
cost_etap = d['cos_theta']

#%%

h_cost_etap = B.histo(cost_etap,    bins = 20,  title = "cost_etap", xlabel = "$cos{\\theta\eta^{'}}$")

h_etappi0 = B.histo(metappi0, range = [1.0, 2.02], bins = 20,  
               xlabel = "$M(\eta'\pi^{0})$", title = "Invariant mass of $\eta^{'}\pi^{0}$ ")

h2_etappi0_GJ = B.histo2d(metappi0, cost_etap, range = [[1.0, 2.02], [-1, 1]],
                bins = (24, 24), title = 'Angular distribution in GJ Frame ', 
                xlabel = "$M(\eta'\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

h2_etappi0_pi0p = B.histo2d(metappi0 , mpi0p,   bins = (30, 30),  title = "etappi0 vs pi0p ", 
                        range = [[1.05, 2.0], [1.05, 2.2]], xlabel = "$M({\eta^{'}\pi^{0}})$", ylabel = "$M({\pi^{0}p})$")
