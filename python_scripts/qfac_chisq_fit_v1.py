#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:09:38 2021

@author: rupeshdotel
"""

#qfactors with least square fitting

import numpy as np
import matplotlib.pyplot as plt
import LT.box as B

#%%

A = B.Parameter(50., 'A')
x0 = B.Parameter(0.956, 'x0')
sigma = B.Parameter(.005, 'sigma')

def gaus(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sigma()**2))


a0 = B.Parameter(1., 'a0')
a1 = B.Parameter(1., 'a1')


def lin_bkg(x):
    return a0() + a1()*x 

def signal(x):
    return gaus(x) + lin_bkg(x)

#%%

f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18_incamo.npz')
metappi0 = f['metappi0']
cost_etap = f['cost_etap']

#%%
# m_etap is the mass of etaprime  
m_etap = f['metap']
# M_min and M_max are the left and right edge values
M_min = 0.90; M_max = 1.02
sel = (M_min < m_etap) & (m_etap < M_max)
# M_inv_sel are the selected events within the range
M_inv_sel = m_etap[sel]

# cos_gj is the costheta in Gottfried Jackson frame
cos_gj = f['cost_etap'][:]

cos_gj_sel = cos_gj[sel]

# cos_gj_range is the total range of costheta variable
cos_gj_range = 2.0

# dcos_gj is the neighborhood 
dcos_gj = 0.05

# Nf are  the number of neighboring events we choose
Nf = 800

# Min number of events per event in neighborhood for warning
Nf_min = 500
# k is an index for testing
k = 1

qf = np.ones_like(M_inv_sel[:])
# X0 are the mean values from the fit
X0 = np.ones_like(M_inv_sel[:])

# S0 are the sigma values from the fit
S0 = np.ones_like(M_inv_sel[:])
do_plot  = True

for i,M_inv_loc in enumerate(M_inv_sel[:k]):
    cos_gj_loc = cos_gj_sel[i]
    # select neghboring events for the current event
    sel_n = np.abs(cos_gj_loc -cos_gj_sel) <= dcos_gj
    
    if sel_n.sum() < Nf_min:
        print(f'found only {sel_n.sum()} events for neighborhood of m = {M_inv_loc:.3f}, theta = {cos_gj_loc:.3f}')
        
        
    M_inv_neighbor = M_inv_sel[sel_n][:Nf]
   
    
    h = B.histo(M_inv_neighbor, bins = 22)
    h_sel = h.bin_content > 0
    M = h.bin_center[h_sel]
    C = h.bin_content[h_sel]
    dC = h.bin_error[h_sel]
    
    A = B.Parameter(C.max(), 'A')
    x0 = B.Parameter(0.956, 'x0')
    sigma = B.Parameter(.005, 'sigma')
    a0 = B.Parameter(1., 'a0')
    a1 = B.Parameter(1., 'a1')
    
    fit = B.genfit(signal, [A, x0, sigma,  a0, a1],
                   x = M, y = C, y_err = dC, print_results = True )
    
    if do_plot:
        plt.figure()
        h.plot_exp()
        B.plot_line(fit.xpl, fit.ypl)
        mr = np.linspace(M[0], M[-1], 1000)
        B.plot_line(mr, gaus(mr))
        B.plot_line(mr, lin_bkg(mr))

    A = B.Parameter(fit.parameters[0].value, 'A')
    x0 = B.Parameter(fit.parameters[1].value, 'x0')
    sigma = B.Parameter(fit.parameters[2].value, 'sigma')
    a0 = B.Parameter(fit.parameters[3].value, 'a0')
    a1 = B.Parameter(fit.parameters[4].value, 'a1')
    #x0 = fit.parameters[1].value
    #sigma = fit.parameters[2].value
    #a0 = fit.parameters[3].value
    #a1 = fit.parameters[4].value



    ws = gaus(M_inv_loc)
    wb = lin_bkg(M_inv_loc)
    # calculate q-factor  
    q = ws/(ws+wb)
    qf[i] = q
    X0[i] = x0.value
    S0[i] = sigma.value

#%%

metap_min = 0.90
metap_max = 1.02
metap_binwidth = 3e-3 # bin width in GeV
bins_metap = int((metap_max - metap_min)/metap_binwidth)
h  = B.histo(M_inv_sel[:], bins = bins_metap)   
hs  = B.histo(M_inv_sel[:], weights = qf,  bins = bins_metap)   
hb  = B.histo(M_inv_sel[:], weights = 1-qf,  bins = bins_metap)   


plt.figure()
h.plot_exp()
hs.plot_exp()
hb.plot_exp()
plt.title("Q FACTOR weighted histogram") 
plt.xlabel("Invariant Mass of etaprime")
plt.ylabel("counts")
    
#%%

cost_min = -1.0
cost_max = 1.0
bins_cost = 20

metappi0_min = 1.0
metappi0_max = 2.5
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max - metappi0_min)/metappi0_binwidth)

h2_ep_cost_weighted = B.histo2d(metappi0[sel][:], cost_etap[sel][:],  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               weights = qf)


#%%

np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qweights.npz',
         qf = qf)




