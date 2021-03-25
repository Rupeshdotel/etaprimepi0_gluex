#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:44:29 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import  matplotlib.pyplot  as plt
#%%

A = B.Parameter(500., 'A')
x0 = B.Parameter(0.956, 'x0')
sig = B.Parameter(.005, 'sigma')

def gaus(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2))


a0 = B.Parameter(1., 'a0')
a1 = B.Parameter(1., 'a1')


def lin_bkg(x):
    return a0() + a1()*x 

def signal(x):
    return gaus(x) + lin_bkg(x)



#%%
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18.npz')

metap = d['metap']
metappi0 = d['metappi0']
cost_etap = d['cost_etap']

#%%

cost_min = -1.0
cost_max = 1.0
bins_cost = 4

metap_min = 0.85
metap_max = 1.05
metap_binwidth = 5e-3 # bin width in GeV
bins_metap = int((metap_max - metap_min)/metap_binwidth)



h2_etap_cost = B.histo2d(metap, cost_etap,  bins = (bins_metap, bins_cost), 
                            range = [[metap_min, metap_max],[cost_min, cost_max]],
                       title = "Angular dependence of eta signal",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$\cos\\theta_{GJ}$" )

metappi0_min = 1.0
metappi0_max = 2.5
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max -metappi0_min)/metappi0_binwidth)   

h2_ep_cost = B.histo2d(metappi0, cost_etap, bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

#%%

def fit_draw_histo(h2d, Fit = False):
    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
    (ax1, ax2), (ax3, ax4) = gs.subplots(sharex='col', sharey='row')
    ax = [ax1, ax2, ax3, ax4]    
    
    for i in range(h2d.nbins_y):
        h_xproj = h2d.project_x(bins=[i])
        ax[i].plot(h_xproj.xpl, h_xproj.ypl)
        if(Fit == True):
            minimum = 0.90
            maximum = 1.04
            sel  = (minimum < h_xproj.bin_center) & (h_xproj.bin_center < maximum)
            
            fit = B.genfit(signal, [A, x0, sig,  a0, a1],
                           x = h_xproj.bin_center[sel] , y = h_xproj.bin_content[sel], y_err = h_xproj.bin_error[sel] )
            ax[i].plot(fit.xpl, fit.ypl)

        
        
        
        
#%%
 

    
    




