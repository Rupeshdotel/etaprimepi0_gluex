#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 20:32:15 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import ROOT as R
from root_numpy import tree2array
import matplotlib.pyplot as plt

#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/accep_mctree.root")

intree = rfile.Get('qfactortree')

d = tree2array(intree)

   
#%%   

mpi0p = d['pi0p']
metappi0 = d['etaprimepi0mass']
cost_etap = d['etaprimecosthetaGJ']

#%%


h_etappi0 = B.histo(metappi0, range = [1.0, 2.02], bins = 20,  
               xlabel = "$M(\eta'\pi^{0})$", title = "Invariant mass of $\eta^{'}\pi^{0}$ ")

h2_etappi0_pi0p = B.histo2d(metappi0 , mpi0p,   bins = (30, 30),  title = "etappi0 vs pi0p ", 
                        range = [[1.05, 2.0], [1.05, 2.2]], xlabel = "$M({\eta^{'}\pi^{0}})$", ylabel = "$M({\pi^{0}p})$")

h2_etappi0_GJ = B.histo2d(metappi0, cost_etap, range = [[1.0, 2.02], [-1, 1]],
                bins = (24, 24), title = 'Angular distribution in GJ Frame ', 
                xlabel = "$M(\eta'\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")