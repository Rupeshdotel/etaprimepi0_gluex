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

f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18.npz')
metappi0 = f['metappi0']
cost_etap = f['cost_etap']

#%%

x = f['metap']
M_min = 0.90; M_max = 1.02
sel = (M_min < x) & (x < M_max)
xsel = x[sel]

yr1 = f['cost_etap'][:]
yr1 = yr1[sel]
y1_range = 2.0
dy = 0.2 
Nf = 800
k = 5000
qf = np.ones_like(xsel[:])
X0 = np.ones_like(xsel[:])
S0 = np.ones_like(xsel[:])


for i,x in enumerate(xsel[:]):
    y11 = yr1[i]
    sel_n = np.sqrt(((y11 - yr1)/y1_range)**2)  <= dy
    xf = xsel[sel_n][:Nf]
   
    h = B.histo(xf, bins = 22)
    M = h.bin_center
    C = h.bin_content
    dC = h.bin_error
    
    A = B.Parameter(C.max(), 'A')
    x0 = B.Parameter(0.956, 'x0')
    sigma = B.Parameter(.005, 'sigma')
    a0 = B.Parameter(1., 'a0')
    a1 = B.Parameter(1., 'a1')
    
    fit = B.genfit(signal, [A, x0, sigma,  a0, a1],
                   x = M, y = C, y_err = dC )
    
    
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



    ws = gaus(x)
    wb = lin_bkg(x)
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
h  = B.histo(xsel[:], bins = bins_metap)   
hs  = B.histo(xsel[:], weights = qf,  bins = bins_metap)   
hb  = B.histo(xsel[:], weights = 1-qf,  bins = bins_metap)   


plt.figure()
h.plot_exp()
hs.plot_exp()
hb.plot_exp()
    
#%%

cost_min = -1.0
cost_max = 1.0
bins_cost = 30

metappi0_min = 1.0
metappi0_max = 2.5
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max - metappi0_min)/metappi0_binwidth)

h2_ep_cost_weighted = B.histo2d(metappi0[sel], cost_etap[sel],  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               weights = qf)


#%%
"""
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/signal_events.npz',
         
         
        qf = qf,  metap = xsel,  metappi0 = metappi0[sel], cost_etap = cost_etap[sel])

"""


