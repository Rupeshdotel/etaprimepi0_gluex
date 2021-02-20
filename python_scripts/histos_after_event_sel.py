#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 13:31:08 2021

@author: rupeshdotel
"""
import matplotlib.pyplot as plt
import numpy as np
import LT.box as B

#%%
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_selected_events_less.npz')

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

t_etap = d['t_etap']
mant = d['mant']
num_unusedshowers = d['num_unusedshowers']

mpipp = d['mpipp']
mpi0p = d['mpi0p']
mpippimpi0 = d['mpippimpi0']

cost_pi0 = d['cost_pi0']
pi0phiGJ = d['pi0phiGJ']

cost_etap = d['cost_etap']
etaprimephiGJ = d['etaprimephiGJ']

pt_p = d['pt_p']
pt_etap = d['pt_etap']
pt_pi0 = d['pt_pi0']

pl_p = d['pl_p']
pl_etap = d['pl_etap']
pl_pi0 = d['pl_pi0']

#%%

h_cost_etap = B.histo(cost_etap,    bins = 20,  title = "cost_etap", xlabel = "$cos{\\theta\eta^{'}}$")
h2_etappi0_pi0p = B.histo2d(metappi0 , mpi0p,   bins = (30, 30),  title = "etappi0 vs pi0p ", 
                        range = [[1.05, 3], [1.05, 3]], xlabel = "$M({\eta^{'}\pi^{0}})$", ylabel = "$M({\pi^{0}p})$")


h_t_etap = B.histo(t_etap,    bins = 20,  title = "t-etaprime", xlabel = "$t_{\eta^{'}}$")
h2_t_etap_mpi0p = B.histo2d(t_etap, mpi0p,   bins = (25, 25),  title = "pi0proton vs t-etaprime ", 
                        range = [[0, 3], [1.05, 3]], xlabel = "$t_{\eta^{'}}$", ylabel = "$M({\pi^{0}p})$")

h_pt_p = B.histo(pt_p,   bins = 20,  title = "transverse mometum of proton in cm frame", xlabel = "$pt_p$")
h_pt_etap = B.histo(pt_etap,   bins = 20,  title = "transverse mometum of $\eta^{'}$ in cm frame", xlabel = "$p{t_{\eta^{'}}}$")
h_pt_pi0 = B.histo(pt_p,   bins = 20,  title = "transverse mometum of $\pi^{0}$ in cm frame", xlabel = "$p{t_{\pi^{0}}}$")

h_pl_p = B.histo(pl_p,   bins = 20,  title = "longitudinal mometum of proton in cm frame", xlabel = "$pl_p$")
h_pl_etap = B.histo(pl_etap,   bins = 20,  title = "longitudinal mometum of $\eta^{'}$ in cm frame", xlabel = "$p{l_{\eta^{'}}}$")
h_pl_pi0 = B.histo(pl_p,   bins = 20,  title = "longitudinal mometum of $\pi^{0}$ in cm frame", xlabel = "$p{l_{\pi^{0}}}$")

hS_GJ = B.histo(metappi0, range = [0.9, 3.1], bins = 20,  
               xlabel = "$M(\eta'\pi^{0})$", title = "Invariant mass of $\eta^{'}\pi^{0}$ ")

h2S_GJ = B.histo2d(metappi0, cost_etap, range = [[0.9, 3.1], [-1, 1]],
                bins = (24, 24), title = 'Angular distribution in GJ Frame ', 
                xlabel = "$M(\eta'\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

h2_etappi0_pippimpi0 = B.histo2d(mpippimpi0, metappi0,  bins = (30, 30) ) #title = "etappi0 vs pippimpi0 ", 
#range = [[1.05, 2], [0.9, 3.1]], xlabel = "$M({\pi^{+}\pi^{-}\pi^{0}})$", ylabel = "$M({\eta^{'}\pi^{0}})$" )


#%%
plt.figure();h_t_etap.plot_exp()
plt.figure();h2_t_etap_mpi0p.plot()
plt.figure();h_pt_p.plot_exp()
plt.figure();h_pt_etap.plot_exp()
plt.figure();h_pt_pi0.plot_exp()

plt.figure();h_pl_p.plot_exp()
plt.figure();h_pl_etap.plot_exp()
plt.figure();h_pl_pi0.plot_exp()







