#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 17:25:31 2021

@author: rupeshdotel
"""

import numpy as np
import ROOT as R
from root_numpy import tree2array
import LT.box as B

#%%
#load accepted MC npz file
mc_ac = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/MC_sel_events.npz')

#load the variables from the file
m_etappi0_ac = mc_ac['metappi0']
cost_etap_ac = mc_ac['cost_etap']


#%%
#load the generated tree and get the variables saved in the tree
rfile = R.TFile('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/gen_mctree.root')
#%%
intree =  rfile.Get('gen_mctree')
d = tree2array(intree)

#%%
m_etappi0_gen = d['metaprimepi0']
cost_etap_gen = d['cos_theta']

#%%
cost_min = -1.0
cost_max = 1.0
bins_cost = 20

metappi0_min = 1.0
metappi0_max = 2.0
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max - metappi0_min)/metappi0_binwidth)

h2_ep_cost_ac = B.histo2d(m_etappi0_ac, cost_etap_ac,  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame")


h2_ep_cost_gen = B.histo2d(m_etappi0_gen, cost_etap_gen,  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame")
#%%

m_gen = B.histo(m_etappi0_gen, bins = 100 )
m_ac = B.histo(m_etappi0_ac, bins = 100 )

m_acceptance = m_ac/m_gen

#%%

cost_etap_gen = B.histo(cost_etap_gen, bins = 100 )
cost_etap_ac = B.histo(cost_etap_ac, bins = 100 )

costheta_acceptance = cost_etap_ac/cost_etap_gen

#%%



acceptance_2d = h2_ep_cost_ac/h2_ep_cost_gen



