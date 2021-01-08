#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:28:26 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B

#%%

#load the npz data file
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_new_eventspromt_17.npz')

#get the variables as numpy arrays
event_num = d['event_num']
kinfit_CL = d['kinfit_CL']
chisq_ndf = d['chisq_ndf']
num_combos = d['num_combos']
combo_number = d['combo_number']
mpi0 = d['mpi0']
meta = d['meta']
metap = d['metap']

metappi0 = d['metappi0']

mpipp = d['mpipp']
mpi0p = d['mpi0p']
mpippimpi0 = d['mpippimpi0']

cost_pi0 = d['pi0costhetaGJ']
pi0phiGJ = d['pi0phiGJ']

cost_etap = d['etaprimecosthetaGJ']
etaprimephiGJ = d['etaprimephiGJ']

#new variables

mpi013 = d['mpi013']
mpi024 = d['mpi024']
mpi014 = d['mpi014']
mpi023 = d['mpi023']

photon1_sq = d['photon1_sq']
photon2_sq = d['photon2_sq']
photon3_sq = d['photon3_sq']
photon4_sq = d['photon4_sq']

mant = d['mant']
num_unusedshowers = d['num_unusedshowers']



#%%
# read the parameters data file 
c_par = B.LT.parameterfile.pfile('parameters2.data')

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

quality_photons = showerqual_1 & showerqual_2 & showerqual_3 & showerqual_4

ex_showers = B.in_between(extrashowers_min, extrashowers_max, num_unusedshowers) # extra showers window

#combine selection windows
sel_win = pi0 & eta & etap & t & quality_photons

#choose 2pi0 veto logic
veto_2pi0 = ~ ( (pi013 & pi024) | (pi014 & pi023) )

#veto omega
veto_omega = ~(omega)

# combine vetos
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


h_t = B.histo(mant[sel],  bins = 40, range = [0, 4],  title = "mometum transfer",
               xlabel = "t")

h_pi0p = B.histo(mpi0p[sel],  bins = 24, title = "$\pi^{0}p$",
               xlabel = "$M(\pi^{0}p)$")

h_pipp = B.histo(mpipp[sel],  bins = 24, title = "$\pi^{+}p$",
               xlabel = "$M(\pi^{+}p$)")

h_us = B.histo(num_unusedshowers[sel],  range = [-0.5, 11.5], bins = 12, title = "Unused showers",
               xlabel = "Number of unused showers")
