#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:21:10 2021

@author: rupeshdotel
"""

""" classical sideband study"""

import matplotlib.pyplot as plt
import numpy as np
import LT.box as B

#%%

#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_17_acc_sub.npz')
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_17.npz')

event_num = d['event_num']
kinfit_CL = d['kinfit_CL']
chisq_ndf = d['chisq_ndf']
num_combos = d['num_combos']
combo_number = d['combo_number']
mpi0 = d['mpi0']
meta = d['meta']
metap = d['metap']

metappi0 = d['metappi0']

mpi013 = d['mpi013']
mpi024 = d['mpi024']
mpi014 = d['mpi014']
mpi023 = d['mpi023']

mant = d['mant']
num_unusedshowers = d['num_unusedshowers']

mpipp = d['mpipp']
mpi0p = d['mpi0p']
mpippimpi0 = d['mpippimpi0']

cost_pi0 = d['cost_pi0']
pi0phiGJ = d['pi0phiGJ']

cost_etap = d['cost_etap']
etaprimephiGJ = d['etaprimephiGJ']

acc_w = d['acc_w']

#%%


mep_bins = 25 # bins etaprime invariant mass

#etaprime histo
h_ep = B.histo(metap,   bins = mep_bins,  title = "$\eta'$", xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

#read the ASCII data file
par = B.LT.parameterfile.pfile('sideband.data')

# load the parameters values
la = par.get_value('la'); lb = par.get_value('lb') # left sideband edges
ra = par.get_value('ra'); rb = par.get_value('rb') # right sideband edges

ma = par.get_value('ma'); mb = par.get_value('mb') # edges for middle region

#select left, middle and right regions from etaprime signal
sel_left = ((la < metap) & (metap < lb))
sel_middle = ((ma < metap) & (metap < mb))
sel_right = ((ra < metap) & (metap < rb))




#%%
#get the bkg counts using sideband technique

Yr = h_ep.sum(ra, rb) # total yield in right sideband
rw = rb - ra #right sideband width

# yield per bin right sideband
yr = Yr[0]/rw

Yl = h_ep.sum(la, lb) # total yield in left sideband
lw = lb - la #left sideband width

# yield per bin left sideband
yl = Yl[0]/lw

"""
#option1 
ys = 0.5 * (yl + yr)
sw = mb - ma #signal width

Yt = h_ep.sum(ma, mb) # total bkg + signal yield in the peak region
Nb = ys * sw  #bkg yield in the peak region
Ys =  Yt - Nb #signal yield in the peak region

"""

#option 2 more general ( better option)
xl = (la + lb)/2
xr = (ra + rb)/2
x = (ma + mb)/2


yb = yl + (yr - yl) * (x - xl) /(xr - xl)
sw = mb - ma #signal width

Yt = h_ep.sum(ma, mb) # total bkg + signal yield in the peak region
Nb = yb * sw  #bkg yield in the peak region
Ys =  Yt - Nb #signal yield in the peak region


#%%

#look at the 2d angular plots in the sidebands and after the bkg is  subtracted using sideband technique

h2_left = B.histo2d(metappi0[sel_left], cost_etap[sel_left],  bins = ( 18, 18 ),
                       title = "Angular distribution in GJ Frame left-sideband",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

h2_right = B.histo2d(metappi0[sel_right], cost_etap[sel_right],  bins = ( 18, 18 ),
                       title = "Angular distribution in GJ Frame right-sideband ",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

h2_middle = B.histo2d(metappi0[sel_middle], cost_etap[sel_middle],  bins = ( 18, 18 ),
                       title = "Angular distribution in GJ Frame middle-band",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )


h2_signal = h2_middle - 0.5 * (h2_left + h2_right)
h2_signal.xlabel =  "$M(\eta^{'}\pi^{0})$"
h2_signal.ylabel = "$\cos\\theta_{GJ}$" 
h2_signal.title  = "Angular distribution in GJ Frame sideband subtracted"
h2_signal.y_bins =  h2_left.y_bins
h2_signal.y_bin_center =  h2_left.y_bin_center


#%%
h2_left = B.histo2d(metappi0[sel_left], cost_etap[sel_left],  bins = ( 18, 18 ), weights = acc_w[sel_left],
                       title = "Angular distribution in GJ Frame left-sideband",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

h2_right = B.histo2d(metappi0[sel_right], cost_etap[sel_right],  bins = ( 18, 18 ), weights = acc_w[sel_right],
                       title = "Angular distribution in GJ Frame right-sideband ",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

h2_middle = B.histo2d(metappi0[sel_middle], cost_etap[sel_middle],  bins = ( 18, 18 ), weights = acc_w[sel_middle],
                       title = "Angular distribution in GJ Frame middle-band",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )


h2_signal = h2_middle - 0.5 * (h2_left + h2_right)
h2_signal.xlabel =  "$M(\eta^{'}\pi^{0})$"
h2_signal.ylabel = "$\cos\\theta_{GJ}$" 
h2_signal.title  = "Angular distribution in GJ Frame sideband subtracted"
h2_signal.y_bins =  h2_left.y_bins
h2_signal.y_bin_center =  h2_left.y_bin_center


#%%

#look at the 1d etaprimepi0 inavriant mass plots in the sidebands and after the bkg is  subtracted using sideband technique

h_left = B.histo(metappi0[sel_left], range = [1.3, 3.0] ,  bins =  15,
                       title = "left-sideband",
               xlabel = "$M(\eta^{'}\pi^{0})$" )

h_middle = B.histo(metappi0[sel_middle],  range = [1.3, 3.0] ,   bins =  15,
                       title = " middle-band",
               xlabel = "$M(\eta^{'}\pi^{0})$" )

h_right = B.histo(metappi0[sel_right],   range = [1.3, 3.0] ,  bins =  15,
                       title = "right-sideband",
               xlabel = "$M(\eta^{'}\pi^{0})$" )


h_signal = h_middle - 0.5 * (h_left + h_right)
h_signal.xlabel =  "$M(\eta^{'}\pi^{0})$ $GeV/c^{2}$"
h_signal.title  = "sideband subtracted Invariant mass distribution"


#%%



h_left = B.histo(metappi0[sel_left],  range = [1.3, 3.0],    bins =  15, weights = acc_w[sel_left],
                       title = "left-sideband",
               xlabel = "$M(\eta^{'}\pi^{0})$" )

h_middle = B.histo(metappi0[sel_middle],  range = [1.3, 3.0],  bins =  15, weights = acc_w[sel_middle],
                       title = " middle-band",
               xlabel = "$M(\eta^{'}\pi^{0})$" )

h_right = B.histo(metappi0[sel_right],  range = [1.3, 3.0],  bins =  15, weights = acc_w[sel_right],
                       title = "right-sideband",
               xlabel = "$M(\eta^{'}\pi^{0})$" )


h_signal = h_middle - 0.5 * (h_left + h_right)
h_signal.xlabel =  "$M(\eta^{'}\pi^{0})$"
h_signal.title  = "sideband subtracted Invariant mass distribution"










