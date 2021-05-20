#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 19:49:36 2021

@author: rupeshdotel

numpy array to root histogram

"""

import numpy as np
import ROOT as R
from root_numpy import fill_hist
import LT.box as B

#%%


d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI.npz')
mpi0 = d['mpi0']
metap = d['metap']
metappi0 = d['metappi0']
cost_etap = d['cost_etap']
cost_eppi0 = np.array([metappi0, cost_etap])
cost_eppi0 = cost_eppi0.T
#combined = np.vstack((metappi0, cost_etap)).T
f_w = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/weights.npz')
qs = f_w['weights']


#%%
R.gROOT.Reset()

h_pi0 = R.TH1F("h_pi0", "$M(\pi^{0})$;" , 80, 0.115,  0.155)
h_ep = R.TH1F("h_ep", "$M(\eta^{\prime})$;" , 80, 0.85,  1.05)


h_ep_pi0 = R.TH1F("h_ep_pi0", "$M(\eta^{\prime}\pi^{0})$;" , 80, 0.9,  3.0)
h_ep_pi0_cost = R.TH2F("h_ep_pi0_cost", "$cos\\theta_{GJ} VS M(\eta^{\prime}\pi^{0})$" , 80, 0.9,  3.0, 80, -1.0, 1.0)

h_ep_pi0_weighted = R.TH1F("h_ep_pi0_weighted", "$M(\eta^{\prime}\pi^{0})$;" , 80, 0.9,  3.0)
h_ep_pi0_cost_weighted = R.TH2F("h_ep_pi0_cost_weighted", "$cos\\theta_{GJ} VS M(\eta^{\prime}\pi^{0})$" , 80, 0.9,  3.0, 80, -1.0, 1.0)

fill_hist(h_pi0, mpi0)
fill_hist(h_ep, metap)

fill_hist(h_ep_pi0, metappi0)
fill_hist(h_ep_pi0_cost, cost_eppi0)
fill_hist(h_ep_pi0_weighted, metappi0, weights = qs)
fill_hist(h_ep_pi0_cost_weighted, cost_eppi0, weights = qs)
outHistFile = R.TFile.Open("/Users/rupeshdotel/analysis/work/pi0pippimeta/for_Malte.root" ,"RECREATE")
outHistFile.cd()

h_pi0.Write()
h_ep.Write()
h_ep_pi0.Write()
h_ep_pi0_cost.Write()

h_ep_pi0_weighted.Write()
h_ep_pi0_cost_weighted.Write()



#%%

metap = d['metap']

#read the ASCII data file
par = B.LT.parameterfile.pfile('/Users/rupeshdotel/analysis/work/gluexgit/data_files/data/sideband.data')

# load the parameters values
la = par.get_value('la'); lb = par.get_value('lb') # left sideband edges
ra = par.get_value('ra'); rb = par.get_value('rb') # right sideband edges

ma = par.get_value('ma'); mb = par.get_value('mb') # edges for middle region

#select left, middle and right regions from etaprime signal
sel_left = ((la < metap) & (metap < lb))
sel_middle = ((ma < metap) & (metap < mb))
sel_right = ((ra < metap) & (metap < rb))

h_ep_pi0_l = R.TH1F("h_ep_pi0_l", "$M(\eta^{\prime}\pi^{0})$;" , 80, 0.9,  3.0)
h_ep_pi0_r = R.TH1F("h_ep_pi0_r", "$M(\eta^{\prime}\pi^{0})$;" , 80, 0.9,  3.0)
h_ep_pi0_m = R.TH1F("h_ep_pi0_m", "$M(\eta^{\prime}\pi^{0})$;" , 80, 0.9,  3.0)

h_ep_pi0_cost_l = R.TH2F("h_ep_pi0_cost_l", "$ cos\\theta_{GJ}  VS M(\eta^{\prime}\pi^{0})$;" , 
                         80, 0.9,  3.0, 80, -1.0, 1.0 )
h_ep_pi0_cost_r = R.TH2F("h_ep_pi0_cost_r", "$ cos\\theta_{GJ}  VS M(\eta^{\prime}\pi^{0})$;" ,
                         80, 0.9,  3.0, 80, -1.0, 1.0)
h_ep_pi0_cost_m = R.TH2F("h_ep_pi0_cost_m", "$ cos\\theta_{GJ}  VS M(\eta^{\prime}\pi^{0})$;" , 
                         80, 0.9,  3.0, 80, -1.0, 1.0)


fill_hist(h_ep_pi0_l, metappi0[sel_left])
fill_hist(h_ep_pi0_r, metappi0[sel_right])
fill_hist(h_ep_pi0_m, metappi0[sel_middle])

cost_eppi0_l = np.array([metappi0[sel_left], cost_etap[sel_left]])
cost_eppi0_r = np.array([metappi0[sel_right], cost_etap[sel_right]])
cost_eppi0_m = np.array([metappi0[sel_middle], cost_etap[sel_middle]])

cost_eppi0_l = cost_eppi0_l.T
cost_eppi0_r = cost_eppi0_r.T
cost_eppi0_m = cost_eppi0_m.T


fill_hist(h_ep_pi0_cost_l, cost_eppi0_l)
fill_hist(h_ep_pi0_cost_r, cost_eppi0_r)
fill_hist(h_ep_pi0_cost_m, cost_eppi0_m)


h_ep_pi0_l.Write()
h_ep_pi0_r.Write()
h_ep_pi0_m.Write()

h_ep_pi0_cost_l.Write()
h_ep_pi0_cost_r.Write()
h_ep_pi0_cost_m.Write()


outHistFile.Close()







