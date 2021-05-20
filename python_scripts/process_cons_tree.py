#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 00:02:31 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import ROOT as R
from root_numpy import tree2array
import matplotlib.pyplot as plt

#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree_for_may_10_gluexI.root")
intree = rfile.Get('qfactortree')
d = tree2array(intree)

#%%



mm2m = d['mm2m']
mpi013 = d['mpi013']
mpi014 = d['mpi014']
mpi023 = d['mpi023']
mpi024 = d['mpi024']


metap = d['metap']
metappi0 = d['metappi0']
cost_etap = d['cos_t']
phi_etap = d['phi_gj']
mpippimpi0 = d['mpippimpi0']
mpi0p = d['mpi0p']
dt = d['dt']
be = d['be']
mant = d['mant']
photon1_sq = d['photon1_sq']
photon2_sq = d['photon2_sq']
photon3_sq = d['photon3_sq']
photon4_sq = d['photon4_sq']

"""
mm2 = d['mm2']
mpi0 = d['mpi0']
meta = d['meta']
mpi0m = d['mpi0m']
metam = d['metam']
metapm = d['metapm']


mpi013m = d['mpi013m']
mpi014m = d['mpi014m']
mpi023m = d['mpi023m']
mpi024m = d['mpi024m']
mpippimpi0m = d['mpippimpi0m']
mpi0pm = d['mpi0pm']
mantm = d['mantm']
event_num = d['event_num']
run_num = d['run_num']

px_pr = d['px_pr']
px_etapr = d['px_etapr']
px_pi0 = d['px_pi0']

py_pr = d['py_pr']
py_etapr = d['py_etapr']
py_pi0 = d['py_pi0']

pz_pr = d['pz_pr']
pz_etapr = d['pz_etapr']
pz_pi0 = d['pz_pi0']

e_pr = d['e_pr']
e_etapr = d['e_etapr']
e_pi0 = d['e_pi0'];

px_beam = d['px_beam']
py_beam = d['py_beam']
pz_beam = d['pz_beam']
e_beam = d['e_beam']

pol = d['pol']
"""

#%%

# read the parameters data file 
c_par = B.LT.parameterfile.pfile('/Users/rupeshdotel/analysis/work/gluexgit/data_files/data/parameters5.data')

#--------get the selection windows from data file-------------#

mm2_min = c_par.get_value('mm2_min'); mm2_max =  c_par.get_value('mm2_max') 
mpi0_min = c_par.get_value('mpi0_min'); mpi0_max =  c_par.get_value('mpi0_max') 
meta_min = c_par.get_value('meta_min'); meta_max =  c_par.get_value('meta_max') 
metap_min = c_par.get_value('metap_min'); metap_max =  c_par.get_value('metap_max') 

veto_2pi0_min = c_par.get_value('twopi0veto_min'); veto_2pi0_max = c_par.get_value('twopi0veto_max') 
mpippimpi0_min = c_par.get_value('mpippimpi0_min'); mpippimpi0_max =  c_par.get_value('mpippimpi0_max')
mpi0p_min = c_par.get_value('mpi0p_min'); mpi0p_max =  c_par.get_value('mpi0p_max') 

dt_min = c_par.get_value('dt_min'); dt_max = c_par.get_value('dt_max')
t_min = c_par.get_value('mant_min'); t_max = c_par.get_value('mant_max')
sq_min = c_par.get_value('sq_min'); sq_max = c_par.get_value('sq_max')



mm2 = B.in_between(mm2_min, mm2_max, mm2m) #missingmass_squared_window
#--------------get the bool arrays for kinfit variables---------#


#pi0 = B.in_between(mpi0_min, mpi0_max,  mpi0) #pi0window
#eta = B.in_between(meta_min, meta_max,  meta) #etawindow
etap = B.in_between(metap_min, metap_max,  metap) #etaprimewindow

pi013 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi013) #2pi0 windows 
pi024 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi024)
pi014 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi014)
pi023 = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi023)

omega = B.in_between(mpippimpi0_min, mpippimpi0_max, mpippimpi0) #omega window 
pi0p = B.in_between(mpi0p_min, mpi0p_max, mpi0p) # pi0p window
t = B.in_between(t_min, t_max,  mant) # momentum transfer window



#-----------time, photon shower quality should be same for kinfit as well as measured variables-----#
showerqual_1 =   B.in_between(sq_min, sq_max,  photon1_sq) # fcal shower quality window
showerqual_2 =   B.in_between(sq_min, sq_max,  photon2_sq)
showerqual_3 =   B.in_between(sq_min, sq_max,  photon3_sq)
showerqual_4 =   B.in_between(sq_min, sq_max,  photon4_sq)

intime = B.in_between(dt_min, dt_max,  dt) # time window

#ex_showers = B.in_between(extrashowers_min, extrashowers_max, num_unusedshowers) # extra showers window
#events that have all photons > 0.5 for shower quality variable


#combine selection windows (bool arrays) in increasing complexity

quality_photons = showerqual_1 & showerqual_2 & showerqual_3 & showerqual_4
#sel_win  = mm2 & intime
#sel_win =  pi0 & eta & etap & intime
#sel_win =  pi0 & eta & etap & intime & t  
#sel_win =  pi0 & eta & etap & intime &  t  
#sel_win =  pi0 & eta & etap &  intime & t & quality_photons 
#sel_win =  pi0 & eta & etap & t & quality_photons & ex_showers
#sel_win =   pi0 & eta & intime & t  
#sel_win =  pi0 & eta & etap &  t & quality_photons 
#sel_win =  intime & t & quality_photons 
#sel_win = pi0 & eta & etap & intime &  t & quality_photons #combine bools for pi0,e eta, etaprime, intime and shower quality
#sel_win = etap & intime &  t & quality_photons
sel_win = quality_photons

veto_2pi0 = ~ ( (pi013 & pi024) | (pi014 & pi023) ) #choose 2pi0 veto logic
veto_omega = ~(omega) #veto omega
#veto_deltap = ~(pi0p) #veto deltap

# combine vetos in increasing complexity
#veto = veto_2pi0
veto = veto_2pi0 & veto_omega
#veto = veto_2pi0 & veto_omega & veto_deltap

# combine veto and selection windows
sel =  sel_win &  veto 
#sel =  veto 
#sel = sel_win


#--------------get the bool arrays for measured variables-------------#
"""
pi0m = B.in_between(mpi0_min, mpi0_max,  mpi0m) #pi0window
etam = B.in_between(meta_min, meta_max,  metam) #etawindow
etapm = B.in_between(metap_min, metap_max,  metapm) #etaprimewindow

pi013m = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi013m) #2pi0 windows 
pi024m = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi024m)
pi014m = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi014m)
pi023m = B.in_between(veto_2pi0_min, veto_2pi0_max, mpi023m)

omegam = B.in_between(mpippimpi0_min, mpippimpi0_max, mpippimpi0m) #omega window 
pi0pm = B.in_between(mpi0p_min, mpi0p_max, mpi0pm) # pi0p window

tm = B.in_between(t_min, t_max,  mantm) # momentum transfer window


#combine selection windows 
#sel_winm = pi0m & etam & etapm & intime & tm & quality_photons 
sel_winm = pi0m & etam & intime & tm & quality_photons 


#get veto bools
veto_2pi0m = ~ ((pi013m & pi024m) | (pi014m & pi023m)) #choose 2pi0 veto logic
veto_omegam = ~(omegam) #veto omega
veto_deltapm = ~(pi0pm) #veto deltap

#combine veto bools in increasing complexity
#veto = veto_2pi0m & veto_omegam
vetom = veto_2pi0m & veto_omegam & veto_deltapm



selm = sel_winm & vetom
"""
#%%

h_mm2m = B.histo(mm2m, bins = 30, title = "missing mass squared" )


h_etap = B.histo(metap[sel], bins = 24, title = "$\eta^{'}$", xlabel = "$M(\pi^{+}\pi^{-}\eta)")

h_etap_pi0 = B.histo2d(metap[sel], mpi0[sel],   bins = 24, title = " pi0 vs etaprime Kinfit variables",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel ="$M(\gamma\gamma)$" )



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



h_sq = B.histo(photon1_sq[sel],  range = [0, 1], bins = 24, title = "Photon1 shower quality",
               xlabel = "Shower quality")


h2_etappi0_pi0p = B.histo2d(metappi0[sel] , mpi0p[sel],   bins = (30, 30),  title = "etappi0 vs pi0p ", 
                        range = [[1.05, 3], [1.05, 3]], xlabel = "$M({\eta^{'}\pi^{0}})$", ylabel = "$M({\pi^{0}p})$")


h_pi0p = B.histo(mpi0p[sel],  bins = 40, title = "$\pi^{0}p$",
               xlabel = "$M(\pi^{0}p)$")

#%%
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_cons.npz',
        px_pr = px_pr[sel],
        px_etapr = px_etapr[sel],
        px_pi0 = px_pi0[sel],

        py_pr = py_pr[sel],
        py_etapr = py_etapr[sel],
        py_pi0 = py_pi0[sel],

        pz_pr = pz_pr[sel],
        pz_etapr = pz_etapr[sel],
        pz_pi0 = pz_pi0[sel],

        e_pr = e_pr[sel],
        e_etapr = e_etapr[sel],
        e_pi0 = e_pi0[sel],
        
        px_beam = px_beam[sel],
        py_beam = py_beam[sel],
        pz_beam = pz_beam[sel],
        e_beam = e_beam[sel],

        pol = pol[sel],
        event_num = event_num[sel],
        run_num = run_num[sel],
        
        mpi0 = mpi0[sel],
        meta = meta[sel],
        metap = metap[sel],

        metappi0 = metappi0[sel],
        cost_etap = cost_etap[sel],
        phi_etap = phi_etap[sel]


)






