#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:28:26 2021

@author: rupeshdotel
"""

"doing event selection before making the events unique makes use of  find duplicates possible in make_events_unique.py "

import numpy as np
import LT.box as B
import ROOT as R
from root_numpy import tree2array
import class_fit as cf
import matplotlib.pyplot as plt
#%%

#load the qtree file
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qtree.root")
intree = rfile.Get('qfactortree')
d = tree2array(intree)

event_num = d['event_num']
kinfit_CL = d['kinfit_CL']
chisq_ndf = d['chisq_ndf']
num_combos = d['num_combos']
combo_number = d['combo_number']
mpi0 = d['pi0mass']
meta = d['etamass']
metap = d['etaprimemass']

metappi0 = d['etaprimepi0mass']

mpi013 = d['pi0mass13']
mpi024 = d['pi0mass24']
mpi014 = d['pi0mass14']
mpi023 = d['pi0mass23']

mant = d['mant']
num_unusedshowers = d['num_unusedshowers']

mpipp = d['pipp']
mpi0p = d['pi0p']
mpippimpi0 = d['pippimpi0']

cost_pi0 = d['pi0costhetaGJ']
pi0phiGJ = d['pi0phiGJ']

cost_etap = d['etaprimecosthetaGJ']
etaprimephiGJ = d['etaprimephiGJ']

photon1_sq = d['photon1_sq']
photon2_sq = d['photon2_sq']
photon3_sq = d['photon3_sq']
photon4_sq = d['photon4_sq']


#%%
# read the parameters data file 
c_par = B.LT.parameterfile.pfile('parameters5.data')

# load the parameters values
mpi0_min = c_par.get_value('mpi0_min'); mpi0_max =  c_par.get_value('mpi0_max') 
meta_min = c_par.get_value('meta_min'); meta_max =  c_par.get_value('meta_max') 
metap_min = c_par.get_value('metap_min'); metap_max =  c_par.get_value('metap_max') 

veto_2pi0_min = c_par.get_value('twopi0veto_min'); veto_2pi0_max = c_par.get_value('twopi0veto_max') 
mpippimpi0_min = c_par.get_value('mpippimpi0_min'); mpippimpi0_max =  c_par.get_value('mpippimpi0_max')

t_min = c_par.get_value('mant_min'); t_max = c_par.get_value('mant_max')

sq_min = c_par.get_value('sq_min'); sq_max = c_par.get_value('sq_max')

extrashowers_min = c_par.get_value('extrashowers_min')
extrashowers_max = c_par.get_value('extrashowers_max')

# select windows for variables

pi0 = B.in_between(mpi0_min, mpi0_max,  mpi0) #pi0window
eta = B.in_between(meta_min, meta_max,  meta) #etawindow
etap = B.in_between(metap_min, metap_max,  metap) #etaprimewindow

pi013 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi013) #2pi0 windows 
pi024 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi024)
pi014 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi014)
pi023 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi023)

omega = B.in_between(mpippimpi0_min, mpippimpi0_max, mpippimpi0) #omega window

t = B.in_between(t_min, t_max,  mant) # momentum transfer window

showerqual_1 =   B.in_between(sq_min, sq_max,  photon1_sq) # fcal shower quality window
showerqual_2 =   B.in_between(sq_min, sq_max,  photon2_sq)
showerqual_3 =   B.in_between(sq_min, sq_max,  photon3_sq)
showerqual_4 =   B.in_between(sq_min, sq_max,  photon4_sq)

ex_showers = B.in_between(extrashowers_min, extrashowers_max, num_unusedshowers) # extra showers window


#events that have all photons > 0.5 for shower quality variable
quality_photons = showerqual_1 & showerqual_2 & showerqual_3 & showerqual_4

#combine selection windows in increasing complexity
#sel_win = pi0 & eta & etap
#sel_win = pi0 & eta & etap & t
sel_win = pi0 & eta & etap & t & quality_photons
#sel_win = pi0 & eta & etap & t & quality_photons & ex_showers

#choose 2pi0 veto logic
veto_2pi0 = ~ ( (pi013 & pi024) | (pi014 & pi023) )

#veto omega
veto_omega = ~(omega)

# combine vetos in increasing complexity
#veto = veto_2pi0
veto = veto_2pi0 & veto_omega

# combine veto and selection windows
sel =  sel_win &  veto 

#%%
#define histograms


h_pi0 = B.histo(mpi0[sel],  bins = 24, title = "$\pi^{0}$",
               xlabel = "$M(\gamma\gamma)$")

h_eta = B.histo(meta[sel],  bins = 24, title = "$\eta$",
               xlabel = "$M(\gamma\gamma)$")

h_etap = B.histo(metap[sel],  bins = 24, title = "$\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

h_etappi0 = B.histo(metappi0[sel],  bins = 24, title = "$\eta'\pi^{0}$",
               xlabel = "$M(\eta^{'}\pi^{0})$")

h2_ep_cost = B.histo2d(metappi0[sel], cost_etap[sel],  bins = ( 24, 24 ),
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )


h2_13_24 = B.histo2d(mpi013, mpi024, range = [[0.08, 0.20], [0.08, 0.20]],  bins = ( 24, 24 ),
                       title = "pi0 pi0 events",
               xlabel = "$M(\gamma_{1}\gamma_{3})$", ylabel = "$M(\gamma_{2}\gamma_{4})$")

h2_14_23 = B.histo2d(mpi014, mpi023, range = [[0.08, 0.20], [0.08, 0.20]],  bins = ( 24, 24 ),
                       title = "pi0 pi0 events",
               xlabel = "$M(\gamma_{1}\gamma_{4})$", ylabel = "$M(\gamma_{2}\gamma_{3})$")

h_omega = B.histo(mpippimpi0[sel],  bins = 40,  range = [0.5, 1.5], title = "$\omega$",
               xlabel = "$M(\pi^{+}\pi^{-}\pi^{0})$")


h2_omega = B.histo2d(mpippimpi0[sel], metap[sel], range = [[0.5, 1.5], [0.9, 1.05]],  bins = ( 40, 40 ),
                       title = "pi0 pi0 events",
               xlabel = "$M(\gamma_{1}\gamma_{3})$", ylabel = "$M(\gamma_{2}\gamma_{4})$")


h_t = B.histo(mant[sel],  bins = 40, range = [0, 4],  title = "mometum transfer",
               xlabel = "t")

h_pi0p = B.histo(mpi0p[sel],  bins = 24, title = "$\pi^{0}p$",
               xlabel = "$M(\pi^{0}p)$")

h_pipp = B.histo(mpipp[sel],  bins = 24, title = "$\pi^{+}p$",
               xlabel = "$M(\pi^{+}p$)")

h_us = B.histo(num_unusedshowers[sel],  range = [-0.5, 11.5], bins = 12, title = "Unused showers",
               xlabel = "Number of unused showers")

h_sq = B.histo(photon1_sq[sel],  range = [0, 1], bins = 24, title = "Photon1 shower quality",
               xlabel = "Shower quality")


#%%

plt.figure();h_pi0.plot_exp()
plt.figure();h_eta.plot_exp()
plt.figure();h_etap.plot_exp()
plt.figure();h_etappi0.plot_exp()
plt.figure();h2_ep_cost.plot()



#%%

"""
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events.npz',
         
         
            event_num = event_num[sel],
            kinfit_CL = kinfit_CL[sel],
            chisq_ndf = chisq_ndf[sel],
            num_combos = num_combos[sel],
            combo_number = combo_number[sel],
            mpi0 = mpi0[sel],
            meta = meta[sel],
            metap = metap[sel],
            
            metappi0 = metappi0[sel],
            
            mpi013 = mpi013[sel],
            mpi024 = mpi024[sel],
            mpi014 = mpi014[sel],
            mpi023 = mpi023[sel],
            
            mant = mant[sel],
            num_unusedshowers = num_unusedshowers[sel],
            
            mpipp = mpipp[sel],
            mpi0p = mpi0p[sel],
            mpippimpi0 = mpippimpi0[sel],
            
            cost_pi0 = cost_pi0[sel],
            pi0phiGJ = pi0phiGJ[sel],
            
            cost_etap = cost_etap[sel],
            etaprimephiGJ = etaprimephiGJ[sel],
            
            photon1_sq = photon1_sq[sel],
            photon2_sq = photon2_sq[sel],
            photon3_sq = photon3_sq[sel],
            photon4_sq = photon4_sq[sel]
            
            )

"""






