#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:35:36 2020

@author: rupeshdotel
"""

import numpy as np
import ROOT as R
from root_numpy import tree2array
#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree.root")
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

pi0costhetaGJ = array['pi0costhetaGJ']
pi0phiGJ = array['pi0phiGJ']

etaprimecosthetaGJ = array['etaprimecosthetaGJ']
etaprimephiGJ = array['etaprimephiGJ']

#%%
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/new_eventsprompt_gluexI.npz', 
         
         
        event_num = event_num,
        kinfit_CL = kinfit_CL,
        chisq_ndf = chisq_ndf,
        num_combos = num_combos,
        combo_number = combo_number,
        mpi0 = mpi0,
        meta = meta,
        metap = metap,
        
        metappi0 = metappi0,
        
        mpi013 = mpi013,
        mpi024 = mpi024,
        mpi014 = mpi014,
        mpi023 = mpi023,
        
        mant = mant,
        num_unusedshowers = num_unusedshowers,
        
        mpipp = mpipp,
        mpi0p = mpi0p,
        mpippimpi0 = mpippimpi0,
        
        pi0costhetaGJ = pi0costhetaGJ,
        pi0phiGJ = pi0phiGJ,
        
        etaprimecosthetaGJ = etaprimecosthetaGJ,
        etaprimephiGJ = etaprimephiGJ 
        )