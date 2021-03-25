#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 11:32:54 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B


#%%
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_cons_17_18_pol_0.npz')
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18.npz')


#%%
metap = d['metap']
metappi0 = d['metappi0']
cost_etap = d['cost_etap']


h_etap = B.histo(metap, bins = 32,  title = "$M(\pi^{+}\pi^{-}\eta)$ Kinfit Variables (Mass constrained trees)")


h_etappi0 = B.histo(metappi0,  bins = 40, 
               title = "$M(\eta^{'}\pi^{0})$ Kinfit Variables  (Mass constrianed trees) ")

h2_ep_cost = B.histo2d(metappi0, cost_etap,  bins = ( 25, 25 ),
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )


#%%
#read the ASCII data file
par = B.LT.parameterfile.pfile('/Users/rupeshdotel/analysis/work/gluexgit/data_files/data/sideband.data')

# load the parameters values
la = par.get_value('la'); lb = par.get_value('lb') # left sideband edges
ma = par.get_value('ma'); mb = par.get_value('mb') # edges for middle region
ra = par.get_value('ra'); rb = par.get_value('rb') # right sideband edges



#left band center
il = (la + lb)/2

#middle band center
im = (ma + mb)/2

#right band center
ir = (ra + rb)/2


#select arrays for left, middle and right
sel_l = (la < metap) & (metap < lb)
sel_m = (ma < metap) & (metap < mb)
sel_r = (ra < metap) & (metap < rb)


bool_Nl = (la < h_etap.bins) & (h_etap.bins < lb)
bool_Nm = (ma < h_etap.bins) & (h_etap.bins < mb)
bool_Nr = (ra < h_etap.bins) & (h_etap.bins < rb)

#numbers of bins in left, middle and right bands
Nl = len(h_etap.bins[bool_Nl])
Nm = len(h_etap.bins[bool_Nm])
Nr = len(h_etap.bins[bool_Nr])

#scale factor
f = (im - il)/(ir - il)




#%%
cost_min = -1.0
cost_max = 1.0
bins_cost = 20

metappi0_min = 1.0
metappi0_max = 2.5
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max - metappi0_min)/metappi0_binwidth)

#define histos in 3 differnet regions
#1D histos

#2D histos
h2_ep_cost_left = B.histo2d(metappi0[sel_l], cost_etap[sel_l],  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

h2_ep_cost_middle = B.histo2d(metappi0[sel_m], cost_etap[sel_m],  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$")

h2_ep_cost_right = B.histo2d(metappi0[sel_r], cost_etap[sel_r],  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

#%%
#Bkg histos
h2_ep_cost_Bkg = f * Nm/Nr * h2_ep_cost_right + Nm/Nl *(1 - f) * h2_ep_cost_left

h2_ep_cost_Bkg.xlabel =  "$M(\eta^{'}\pi^{0})$"
h2_ep_cost_Bkg.ylabel = "$\cos\\theta_{GJ}$" 
h2_ep_cost_Bkg.title  = "Bkg histo"
h2_ep_cost_Bkg.nbins_y = h2_ep_cost_left.nbins_y
h2_ep_cost_Bkg.y_bins =  h2_ep_cost_left.y_bins
h2_ep_cost_Bkg.y_bin_center =  h2_ep_cost_left.y_bin_center

#correctedd histos
h2_ep_cost_corr = h2_ep_cost_middle - h2_ep_cost_Bkg



h2_ep_cost_corr.xlabel =  "$M(\eta^{'}\pi^{0})$"
h2_ep_cost_corr.ylabel = "$\cos\\theta_{GJ}$" 
h2_ep_cost_corr.title  = "Angular distribution in GJ Frame sideband subtracted"
h2_ep_cost_corr.nbins_y = h2_ep_cost_left.nbins_y
h2_ep_cost_corr.y_bins =  h2_ep_cost_left.y_bins
h2_ep_cost_corr.y_bin_center =  h2_ep_cost_left.y_bin_center




