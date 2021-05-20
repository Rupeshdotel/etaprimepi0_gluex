#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:00:11 2021

@author: rupeshdotel
"""

import numpy as np
import matplotlib.pyplot as plt
import LT.box as B
import class_fit as cf

#%%
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_17.npz')

event_num = d['event_num']
kinfit_CL = d['kinfit_CL']
chisq_ndf = d['chisq_ndf']
num_combos = d['num_combos']
combo_number = d['combo_number']
mpi0 = d['mpi0']
meta = d['meta']
metap = d['metap']

metappi0 = d['metappi0']

mpi013 = d['mpi013']
mpi024 = d['mpi024']
mpi014 = d['mpi014']
mpi023 = d['mpi023']

mant = d['mant']
num_unusedshowers = d['num_unusedshowers']

mpipp = d['mpipp']
mpi0p = d['mpi0p']
mpippimpi0 = d['mpippimpi0']

cost_pi0 = d['cost_pi0']
pi0phiGJ = d['pi0phiGJ']

cost_etap = d['cost_etap']
etaprimephiGJ = d['etaprimephiGJ']
#%%

mep_bins = 18
mep_min = 0.86
mep_max = 1.05
h = B.histo( metap,  bins = mep_bins,
                   range = [mep_min, mep_max] , title = 'etaprime', xlabel = "$M(\pi^{+}\pi^{-}\eta)$",
                   )
                                  
M = h.bin_center
C = h.bin_content
dC = h.bin_error
mr = np.linspace(M[0], M[-1], 1000)
f = cf.gauss_fit()

# set initial parameters for fitting
f.A.set(C.max()) 
f.x0.set(0.956)
f.sigma.set(0.01)
f.b0.set(0.05)
f.db0.set(0.5)
f.c0.set(500)

# set bounds for the fit parameters
f.A_min.set(0.); f.A_max.set(1e5)
f.x0_min.set(0.95); f.x0_max.set(0.96)
f.sigma_min.set(0.008); f.sigma_max.set(0.016)
f.c0_min.set(0.); f.c0_max.set(1e5)
f.b0_min.set(0.00); f.b0_max.set(0.20)


f.set_fit_list(fit = ['A', 'x0', 'sigma',  'c0', 'b0'])
f.fit_gaussbt(M, C, dC) # gauss peak with bernstein bkg

h.plot_exp()
f.plot_fit()
B.plot_line(mr, f.gauss(mr))
B.plot_line(mr, f.bt_bkg(mr))
             
#%%
s = np.abs(f.gauss(metap))
b = np.abs(f.bt_bkg(metap))
qs_1d = s/(s+b)


hs = B.histo( metap,  bins = mep_bins,
                   range = [mep_min, mep_max] , title = 'etaprime', xlabel = "$M(\pi^{+}\pi^{-}\eta)$",
                 weights = qs_1d  )

hb = B.histo( metap,  bins = mep_bins,
                   range = [mep_min, mep_max] , title = 'etaprime', xlabel = "$M(\pi^{+}\pi^{-}\eta)$",
                 weights = 1 - qs_1d  )      


hS_GJ = B.histo(metappi0, range = [1.3, 3.0], bins = 15,  weights = qs_1d, 
               xlabel = "$M(\eta'\pi^{0})$", title = "Invariant mass of $\eta^{'}\pi^{0}$ Signal Response")



                    
                                   