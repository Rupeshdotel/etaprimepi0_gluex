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

#%%


d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI.npz')
metappi0 = d['metappi0']
cost_etap = d['cost_etap']
cost_eppi0 = np.array([metappi0, cost_etap])
cost_eppi0 = cost_eppi0.T
#combined = np.vstack((metappi0, cost_etap)).T

#%%
R.gROOT.Reset()
h_ep_pi0 = R.TH1F(" h_ep_pi0 ", "$M(\eta^{\prime}\pi^{0})$;" , 80, 0.9,  3.0)
h_ep_pi0_cost = R.TH2F(" h_ep_pi0_cost", "$cos\\theta_{GJ} VS M(\eta^{\prime}\pi^{0})$" , 80, 0.9,  3.0, 80, -1.0, 1.0)

fill_hist(h_ep_pi0, metappi0)
fill_hist(h_ep_pi0_cost, cost_eppi0)
outHistFile = R.TFile.Open("hist.root" ,"RECREATE")
outHistFile.cd()
h_ep_pi0.Write()
h_ep_pi0_cost.Write()
outHistFile.Close()