#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 15:23:52 2020

@author: rupeshdotel
"""

import ROOT as R
from root_numpy import tree2array
#import matplotlib.pyplot as plt
import numpy as np

#%%
#get the input root file with tree
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree_cutexp.root")
intree = rfile.Get('qfactortree')

#%%
array = tree2array(intree)


# get the variables for q_factor analysis

event_num = array['event_num']
kinfit_CL = array['kinfit_CL']
chisq_ndf = array['chisq_ndf']
num_combos = array['num_combos']
combo_number = array['combo_number']
mpi0 = array['pi0mass']
meta = array['etamass']
metap = array['etaprimemass']

metappi0 = array['etaprimepi0mass']

mpi013 = array['pi0mass13']
mpi024 = array['pi0mass24']
mpi014 = array['pi0mass14']
mpi023 = array['pi0mass23']

mant = array['mant']
num_unusedshowers = array['num_unusedshowers']

mpipp = array['pipp']
mpi0p = array['pi0p']
mpippimpi0 = array['pippimpi0']


#%%
#2pi0 veto
mpi013_min = 0.11; mpi013_max = 0.17
mpi024_min = 0.11; mpi024_max = 0.17
mpi014_min = 0.11; mpi014_max = 0.17
mpi023_min = 0.11; mpi023_max = 0.17

#pi0window
mpi0_min = 0.12; mpi0_max = 0.15

#etawindow
meta_min = 0.48; meta_max = 0.60


#etaprimewindow
metap_min = 0.86; metap_max = 1.04


#possible baryons
mpipp_min = 1.35; 
mpi0p_min = 1.35; 

#omega window to be rejected
momega_min = 0.75; momega_max = 0.85; 

#maximum allowed  extra showers
n_maxextrashowers = 1

#maximum mandlestam t allowed
mant_max = 1.




cuts = ((mpi0_min <  mpi0) &  (mpi0  < mpi0_max)) \
        & ((meta_min < meta) & (meta < meta_max)) \
         &   (~((mpi013_min <  mpi013) &  (mpi013  < mpi013_max))) \
         &   (~((mpi024_min <  mpi024) &  (mpi024  < mpi024_max))) \
         &   (~((mpi014_min <  mpi014) &  (mpi014  < mpi014_max))) \
         &    (~((mpi023_min <  mpi023) &  (mpi023  < mpi024_max))) \
         &    (~((momega_min <  mpippimpi0) &  (mpippimpi0  < momega_max))) \
         &    (mpipp > mpipp_min ) \
         &    (mpi0p > mpi0p_min ) \
         &    (num_unusedshowers < n_maxextrashowers )\
         &    (mant  < mant_max) \
         & ((metap_min < metap) & (metap < metap_max)) 






