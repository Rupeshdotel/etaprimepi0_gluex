#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 15:23:52 2020

@author: rupeshdotel
"""

#import ROOT as R
#from root_numpy import tree2array
#import matplotlib.pyplot as plt
import numpy as np
import LT.box as B
import class_fit as cf
from scipy import integrate
#%%
#get the input root file with tree
#rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree_cutexp.root")
#intree = rfile.Get('qfactortree')


d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_new_eventspromt_gluexI.npz')

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

cost_pi0 = d['pi0costhetaGJ']
pi0phiGJ = d['pi0phiGJ']

cost_etap = d['etaprimecosthetaGJ']
etaprimephiGJ = d['etaprimephiGJ']



#%%
#2pi0 veto
twopi0veto_min = 0.11; twopi0veto_max = 0.17
mpi013_min = twopi0veto_min; mpi013_max = twopi0veto_max
mpi024_min = twopi0veto_min; mpi024_max = twopi0veto_max
mpi014_min = twopi0veto_min; mpi014_max = twopi0veto_max
mpi023_min = twopi0veto_min; mpi023_max = twopi0veto_max

#pi0window
mpi0_min = 0.12; mpi0_max = 0.15

#etawindow
meta_min = 0.48; meta_max = 0.60


#etaprimewindow
metap_min = 0.86; metap_max = 1.05


#possible baryons
mpipp_min = 1.35; 
mpi0p_min = 1.35; 

#omega window to be rejected
momega_min = 0.75; momega_max = 0.85; 

#maximum allowed  extra showers
n_maxextrashowers = 3 #optimum value from FOM is 3

#maximum mandlestam t allowed
mant_max = 10.0 # not sure if I need to cut on this 




cuts = ((mpi0_min <  mpi0) &  (mpi0  < mpi0_max)) \
        & ((meta_min < meta) & (meta < meta_max)) \
        & ((metap_min < metap) & (metap < metap_max))\
        &   (~((mpi013_min <  mpi013) &  (mpi013  < mpi013_max))) \
         &   (~((mpi024_min <  mpi024) &  (mpi024  < mpi024_max))) \
         &   (~((mpi014_min <  mpi014) &  (mpi014  < mpi014_max))) \
         &    (~((mpi023_min <  mpi023) &  (mpi023  < mpi024_max))) \
             &    (~((momega_min <  mpippimpi0) &  (mpippimpi0  < momega_max))) \
       &    (num_unusedshowers < n_maxextrashowers )\
           &    (mpi0p > mpi0p_min ) # this cut gives better fit 
        #&   (mpipp > mpipp_min ) # not a very good cut
        
         #&    (mant  < mant_max)\
        


#%%
h2 = B.histo2d(metappi0[cuts], cost_etap[cuts],
               range = [[1.0, 3.0], [-1., 1]] ,bins = [30, 30], 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$", 
                title = 'Angular Distribution in Gottfried-Jackson Frame ')


h_etap_pi0 = B.histo(metappi0[cuts], range = [1.0, 3.0], bins = 40, 
            xlabel = "$M(\eta^{'}\pi^{0})$", title = "Invariant Mass" )

h_ep = B.histo(metap[cuts], range = [0.87, 1.04], bins = 20, 
            xlabel = "$M(\pi^{+}\pi^{-}\eta)$", title = "Invariant Mass" )
#%%
h = h_ep
M = h.bin_center
C = h.bin_content
dC = h.bin_error
mr = np.linspace(M[0], M[-1], 1000)
f = cf.gauss_fit()

# set initial parameters for fitting
f.A.set(C.max()) 
f.x0.set(0.956)
f.sigma.set(0.01)
f.b0.set(0.05)
f.db0.set(0.5)
f.c0.set(5000)


# set bounds for the fit parameters
f.A_min.set(0.); f.A_max.set(1e5)
f.x0_min.set(0.94); f.x0_max.set(0.98)
f.sigma_min.set(0.008); f.sigma_max.set(0.016)
f.c0_min.set(0.); f.c0_max.set(1e5)
f.b0_min.set(0.00); f.b0_max.set(0.10)


f.set_fit_list(fit = ['A', 'x0', 'sigma',  'c0', 'b0'])
f.fit_gaussbt(M, C, dC) # gauss peak with bernstein bkg
h_ep.plot_exp()
B.plot_line(mr, f.signal_bt_bkg(mr))
B.plot_line(mr, f.bt_bkg(mr))
B.plot_line(mr, f.gauss(mr))


l = 0.87
u = 1.04
dM = M[1] - M[0] # get bin width
T_fit = integrate.quad(f.signal_bt_bkg, l, u)[0]/dM
B_fit = integrate.quad(f.bt_bkg, l, u)[0]/dM
S_fit = integrate.quad(f.gauss, l, u)[0]/dM
print(f"\nTotal fitted counts = {T_fit:.0f},  Signal Counts = {S_fit:.0f}, Back Ground = {B_fit:.0f} ")

# total exp. counts

C_tot = C.sum()

print(f"\nExp. to Fit ratio = {C_tot/T_fit:.3f}")

def FOM(S,B):
    return S/(np.sqrt(S+B))

print(f"\nfigure of merit = {FOM(S_fit,B_fit):.4f}")