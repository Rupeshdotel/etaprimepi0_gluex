#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 20:37:41 2021

@author: rupeshdotel
"""

import numpy as np

#%%


fn = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_qfactor.npz')

#%%

px_pr = fn['px_pr']
px_etapr = fn['px_etapr']
px_pi0 = fn['px_pi0']

py_pr = fn['py_pr']
py_etapr = fn['py_etapr']
py_pi0 = fn['py_pi0']

pz_pr = fn['pz_pr']
pz_etapr = fn['pz_etapr']
pz_pi0 = fn['pz_pi0']

e_pr = fn['e_pr']
e_etapr = fn['e_etapr']
e_pi0 = fn['e_pi0']

px_beam = fn['px_beam']
py_beam = fn['py_beam']
pz_beam = fn['pz_beam']
e_beam = fn['e_beam']

pol = fn['pol']
qf = fn['qf']

#m_etap = f['metap']

#M_min = 0.90; M_max = 1.02
#sel = (M_min < m_etap) & (m_etap < M_max)


#%%
"""
px_pr = px_pr[sel]
px_etapr = px_etapr[sel]
px_pi0 = px_pi0[sel]

py_pr = py_pr[sel]
py_etapr = py_etapr[sel]
py_pi0 = py_pi0[sel]

pz_pr = pz_pr[sel]
pz_etapr = pz_etapr[sel]
pz_pi0 = pz_pi0[sel]

e_pr = e_pr[sel]
e_etapr = e_etapr[sel]
e_pi0 = e_pi0[sel]

px_beam = px_beam[sel]
py_beam = py_beam[sel]
pz_beam = pz_beam[sel]
e_beam = e_beam[sel]

pol = pol[sel]
"""

#%%
#qfile = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qweights_18_08_v2.npz')
#qf = qfile['qf']

#%%
sel_amo = pol == -1
sel_0 = pol == 0
sel_90 = pol == 90
sel_45 = pol == 45
sel_135 = pol == 135


#i is the polarization angle, can be 0, 45, 90 135
#i = sel_0 | sel_90 | sel_45 | sel_135
i = sel_amo
pol = pol[i]

px_pr = px_pr[i]
px_etapr = px_etapr[i]
px_pi0 = px_pi0[i]

py_pr = py_pr[i]
py_etapr = py_etapr[i]
py_pi0 = py_pi0[i]

pz_pr = pz_pr[i]
pz_etapr = pz_etapr[i]
pz_pi0 = pz_pi0[i]

e_pr = e_pr[i]
e_etapr = e_etapr[i]
e_pi0 = e_pi0[i]

px_beam = px_beam[i]
py_beam = py_beam[i]
pz_beam = pz_beam[i]
e_beam = e_beam[i]
qf = qf[i]

#%%

np.savez("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_pol_amo.npz",
         pol = pol,

        px_pr = px_pr,
        px_etapr = px_etapr,
        px_pi0 = px_pi0,
        
        py_pr = py_pr,
        py_etapr = py_etapr,
        py_pi0 = py_pi0,
        
        pz_pr = pz_pr,
        pz_etapr = pz_etapr,
        pz_pi0 = pz_pi0,
        
        e_pr = e_pr,
        e_etapr = e_etapr,
        e_pi0 = e_pi0,
        
        px_beam = px_beam,
        py_beam = py_beam,
        pz_beam = pz_beam,
        e_beam = e_beam,
        qf = qf
         
         )










