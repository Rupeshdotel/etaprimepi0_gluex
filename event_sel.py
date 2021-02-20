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
import matplotlib.pyplot as plt
#%%

#load the qtree file
#rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qtree.root")
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree_17.root")
#rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/accep_mctree.root")
#rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/gen_mctree.root")
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

acc_w = d['time_weights']

cost_pi0_GJpi0p = d['cost_pi0']

dt = d['dt']
mis_mass2_m = d['mis_mass2_m']

'''
mis_mass2_m = d['mis_mass2_m']
t_etap = d['t_etap']
pt_p = d['pt_p']; 
pt_etap = d['pt_etap'];
pt_pi0 = d['pt_pi0'];  

pl_p = d['pl_p']; 
pl_etap = d['pl_etap'];
pl_pi0 = d['pl_pi0'];  
'''

#%%
# read the parameters data file 
c_par = B.LT.parameterfile.pfile('parameters5.data')

# load the parameters values
mm2_min = c_par.get_value('mm2_min'); mm2_max =  c_par.get_value('mm2_max') 

mpi0_min = c_par.get_value('mpi0_min'); mpi0_max =  c_par.get_value('mpi0_max') 
meta_min = c_par.get_value('meta_min'); meta_max =  c_par.get_value('meta_max') 
metap_min = c_par.get_value('metap_min'); metap_max =  c_par.get_value('metap_max') 

veto_2pi0_min = c_par.get_value('twopi0veto_min'); veto_2pi0_max = c_par.get_value('twopi0veto_max') 
mpippimpi0_min = c_par.get_value('mpippimpi0_min'); mpippimpi0_max =  c_par.get_value('mpippimpi0_max')
meta_for_omega_min = c_par.get_value('meta_for_omega_min'); meta_for_omega_max = c_par.get_value('meta_for_omega_max') 

dt_min = c_par.get_value('dt_min'); dt_max = c_par.get_value('dt_max')


t_min = c_par.get_value('mant_min'); t_max = c_par.get_value('mant_max')

sq_min = c_par.get_value('sq_min'); sq_max = c_par.get_value('sq_max')

mpi0p_min = c_par.get_value('mpi0p_min'); mpi0p_max =  c_par.get_value('mpi0p_max') 

extrashowers_min = c_par.get_value('extrashowers_min')
extrashowers_max = c_par.get_value('extrashowers_max')

# select windows for variables

#mm2 = B.in_between(mm2_min, mm2_max, mis_mass2_m) #missingmass_squared_window
pi0 = B.in_between(mpi0_min, mpi0_max,  mpi0) #pi0window
eta = B.in_between(meta_min, meta_max,  meta) #etawindow
etap = B.in_between(metap_min, metap_max,  metap) #etaprimewindow

pi013 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi013) #2pi0 windows 
pi024 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi024)
pi014 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi014)
pi023 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi023)

omega = B.in_between(mpippimpi0_min, mpippimpi0_max, mpippimpi0) #omega window
#eta_for_omega = B.in_between(meta_for_omega_min, meta_for_omega_max,  meta) #etawindow

t = B.in_between(t_min, t_max,  mant) # momentum transfer window

intime = B.in_between(dt_min, dt_max,  dt) # time window

showerqual_1 =   B.in_between(sq_min, sq_max,  photon1_sq) # fcal shower quality window
showerqual_2 =   B.in_between(sq_min, sq_max,  photon2_sq)
showerqual_3 =   B.in_between(sq_min, sq_max,  photon3_sq)
showerqual_4 =   B.in_between(sq_min, sq_max,  photon4_sq)

ex_showers = B.in_between(extrashowers_min, extrashowers_max, num_unusedshowers) # extra showers window


#events that have all photons > 0.5 for shower quality variable
quality_photons = showerqual_1 & showerqual_2 & showerqual_3 & showerqual_4

pi0p = B.in_between(mpi0p_min, mpi0p_max, mpi0p) # pi0p window


#combine selection windows in increasing complexity
#sel_win =  pi0 & eta & etap & intime
#sel_win =  pi0 & eta & etap & intime & t  
sel_win =  pi0 & eta & etap & intime &  t  
#sel_win =  pi0 & eta & etap &  intime & t & quality_photons 
#sel_win =  pi0 & eta & etap & t & quality_photons & ex_showers
#sel_win =   pi0 & eta & intime & t  
#sel_win =  pi0 & eta & etap &  t & quality_photons 

#choose 2pi0 veto logic
veto_2pi0 = ~ ( (pi013 & pi024) | (pi014 & pi023) )

#veto omega

veto_omega = ~(omega)

#veto deltap
veto_deltap = ~(pi0p)

# combine vetos in increasing complexity
veto = veto_2pi0
#veto = veto_2pi0 & veto_omega
#veto = veto_2pi0 & veto_omega & veto_deltap

# combine veto and selection windows
#sel  = sel_win
sel =  sel_win &  veto 
#sel =  veto 

#%%
#define histograms

h_pi0 = B.histo(mpi0[sel],  bins = 24, title = "$\pi^{0}$",
               xlabel = "$M(\gamma\gamma)$")

h_eta = B.histo(meta[sel],  bins = 24, title = "$\eta$",
               xlabel = "$M(\gamma\gamma)$")

h_etap = B.histo(metap[sel],  bins = 24, title = "$\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

h_etappi0 = B.histo(metappi0[sel],  bins = 40, title = "$\eta'\pi^{0}$",
               xlabel = "$M(\eta^{'}\pi^{0})$")

h2_ep_cost = B.histo2d(metappi0[sel], cost_etap[sel],  bins = ( 25, 25 ),
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )


h2_13_24 = B.histo2d(mpi013, mpi024, range = [[0.00, 1.0], [0.00, 1.0]],  bins = ( 40, 40 ),
                       title = "pi0 pi0 events",
               xlabel = "$M(\gamma_{1}\gamma_{3})$", ylabel = "$M(\gamma_{2}\gamma_{4})$")

h2_14_23 = B.histo2d(mpi014, mpi023, range = [[0.08, 0.20], [0.08, 0.20]],  bins = ( 24, 24 ),
                       title = "pi0 pi0 events",
               xlabel = "$M(\gamma_{1}\gamma_{4})$", ylabel = "$M(\gamma_{2}\gamma_{3})$")

h_omega = B.histo(mpippimpi0[sel],  bins = 40,  range = [0.5, 1.5], title = "$\omega$",
               xlabel = "$M(\pi^{+}\pi^{-}\pi^{0})$")


h2_omega_metap = B.histo2d(mpippimpi0[sel], metap[sel], range = [[0.5, 1.5], [0.9, 1.05]],  bins = ( 40, 40 ),
                       title = "pi0 pi0 events",
               xlabel = "$M(\gamma_{1}\gamma_{3})$", ylabel = "$M(\gamma_{2}\gamma_{4})$")

h2_omega_metappi0 = B.histo2d(mpippimpi0[sel], metappi0[sel], range = [[0.5, 1.5], [0.9, 2.05]],  bins = ( 40, 40 ),
                       title = "etaprimepi0 vs pippimpi0 events",
               xlabel = "$M(\pi^{+}\pi^{-}\pi^{0})$", ylabel =  "$M(\eta^{'}\pi^{0})$")


h2_omega_meta = B.histo2d(mpippimpi0[sel], meta[sel], range = [[0.5, 1.5], [0.0, 1.05]],  bins = ( 40, 40 ),
                       title = " eta vs pippimpi0 events",
               xlabel = "$M(\pi^{+}\pi^{-}\pi^{0})$", ylabel = "$M(\gamma_{3}\gamma_{4})$")



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


h2_etappi0_pi0p = B.histo2d(metappi0[sel] , mpi0p[sel],   bins = (30, 30),  title = "etappi0 vs pi0p ", 
                        range = [[1.05, 3], [1.05, 3]], xlabel = "$M({\eta^{'}\pi^{0}})$", ylabel = "$M({\pi^{0}p})$")


h_pi0p = B.histo(mpi0p[sel],  bins = 40, title = "$\pi^{0}p$",
               xlabel = "$M(\pi^{0}p)$")

h2_pi0p_cost = B.histo2d(mpi0p[sel], cost_pi0_GJpi0p[sel],  bins = ( 25, 25 ),
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\pi^{0}p)$", ylabel = "$\cos\\theta_{GJ}$" )


#%%

plt.figure();h_pi0.plot_exp()
plt.figure();h_eta.plot_exp()
plt.figure();h_etap.plot_exp()
plt.figure();h_etappi0.plot_exp()
plt.figure();h2_ep_cost.plot()



#%%


np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_17_acc_sub.npz',
         
         
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
            photon4_sq = photon4_sq[sel],
            
           acc_w = acc_w[sel]
            
            
            )

'''
t_etap = t_etap[sel]
pt_p = pt_p[sel]
pt_etap = pt_etap[sel]
pt_pi0 = pt_pi0[sel]

pl_p = pl_p[sel]
pl_etap = pl_etap[sel]
pl_pi0 = pl_pi0[sel]
'''





