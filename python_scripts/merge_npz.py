#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 18:32:04 2021

@author: rupeshdotel
"""

import numpy as np

#%%

"""  

merge S17, S18, and Fall18 

example to concatenate different arrays of different files to a single array and store in a new file
a = np.array([1,2])
b = np.array([3,4])

np.savez('blah.npz', a = a)
np.savez('bla.npz', b = b)

f1 = np.load('blah.npz')
f2 = np.load('bla.npz')

c = np.concatenate([f1['a'], f2['b']])
np.savez('blahbla.npz', c = c )

"""


f17 = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_17_pwa.npz')
f18 = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_18_pwa.npz')
f18_08 = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_18_08_pwa.npz')

#%%

mpi0 = np.concatenate([f17['mpi0'], f18['mpi0'],f18_08['mpi0']])
meta = np.concatenate([f17['meta'], f18['meta'],f18_08['meta']])
metap = np.concatenate([f17['metap'], f18['metap'],f18_08['metap']])
metappi0 = np.concatenate([f17['metappi0'], f18['metappi0'],f18_08['metappi0']])
cost_etap = np.concatenate([f17['cost_etap'], f18['cost_etap'],f18_08['cost_etap']])
#dt = np.concatenate([f17['dt'], f18['dt'],f18_08['dt']])

px_pr = np.concatenate([f17['px_pr'], f18['px_pr'],f18_08['px_pr']])
px_etapr = np.concatenate([f17['px_etapr'], f18['px_etapr'],f18_08['px_etapr']])
px_pi0 = np.concatenate([f17['px_pi0'], f18['px_pi0'],f18_08['px_pi0']])

py_pr = np.concatenate([f17['py_pr'], f18['py_pr'],f18_08['py_pr']])
py_etapr = np.concatenate([f17['py_etapr'], f18['py_etapr'],f18_08['py_etapr']])
py_pi0 = np.concatenate([f17['py_pi0'], f18['py_pi0'],f18_08['py_pi0']])


pz_pr = np.concatenate([f17['pz_pr'], f18['pz_pr'],f18_08['pz_pr']])
pz_etapr = np.concatenate([f17['pz_etapr'], f18['pz_etapr'],f18_08['pz_etapr']])
pz_pi0 = np.concatenate([f17['pz_pi0'], f18['pz_pi0'],f18_08['pz_pi0']])


e_pr = np.concatenate([f17['e_pr'], f18['e_pr'],f18_08['e_pr']])
e_etapr = np.concatenate([f17['e_etapr'], f18['e_etapr'],f18_08['e_etapr']])
e_pi0 = np.concatenate([f17['e_pi0'], f18['e_pi0'],f18_08['e_pi0']])

px_beam = np.concatenate([f17['px_beam'], f18['px_beam'],f18_08['px_beam']])
py_beam = np.concatenate([f17['py_beam'], f18['py_beam'],f18_08['py_beam']])
pz_beam = np.concatenate([f17['pz_beam'], f18['pz_beam'],f18_08['pz_beam']])
e_beam = np.concatenate([f17['e_beam'], f18['e_beam'],f18_08['e_beam']])

#%%





np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI.npz', 
                       mpi0 = mpi0,
                       meta = meta,
                       metap = metap,
                       metappi0 = metappi0,
                       cost_etap = cost_etap,
                       
                       px_pr = px_pr,
                       px_etapr = px_etapr,
                       px_pi0  = px_pi0 ,
                       
                       py_pr = py_pr,
                       py_etapr = py_etapr,
                       py_pi0  = py_pi0 ,
                       
                       pz_pr = pz_pr,
                       pz_etapr = pz_etapr,
                       pz_pi0  = pz_pi0 ,
         
                       e_pr = e_pr,
                       e_etapr = e_etapr,
                       e_pi0  = e_pi0 ,
                       
                      
           
                       px_beam = px_beam,
                       py_beam = py_beam,
                       pz_beam = pz_beam,
                       e_beam = e_beam
           
         
           )


