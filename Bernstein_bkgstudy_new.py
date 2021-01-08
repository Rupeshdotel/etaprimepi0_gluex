#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 11:20:09 2021

@author: rupeshdotel
"""

import matplotlib.pyplot as plt
import numpy as np
import LT.box as B
import class_fit as cf

#%%

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
n_maxextrashowers = 3 #optimum value from FOM is 3 (excludes upperbound)

#maximum mandlestam t allowed
mant_max = 10.0 # not sure if I need to cut on this 



#make a selection of events
sel = ((mpi0_min <  mpi0) &  (mpi0  < mpi0_max)) \
        & ((meta_min < meta) & (meta < meta_max)) \
        & ((metap_min < metap) & (metap < metap_max))\
        &   (~((mpi013_min <  mpi013) &  (mpi013  < mpi013_max))) \
         &   (~((mpi024_min <  mpi024) &  (mpi024  < mpi024_max))) \
         &   (~((mpi014_min <  mpi014) &  (mpi014  < mpi014_max))) \
         &    (~((mpi023_min <  mpi023) &  (mpi023  < mpi024_max))) \
             &    (~((momega_min <  mpippimpi0) &  (mpippimpi0  < momega_max))) \
       &    (num_unusedshowers < n_maxextrashowers )\
           &    (mant  < mant_max)\
               &    (mpi0p > mpi0p_min ) \
        #&   (mpipp > mpipp_min ) 
   
#%%

mep_bins =  20 # mass etaprime bins
ct_bins = 20 # costheta GJ  bins

h2_metap_ct = B.histo2d( metap[sel], cost_etap[sel], bins = (mep_bins, ct_bins),
                   range=np.array([ (0.86, 1.04), (-1.0, 1.0)]), title = '2D')

epm = h2_metap_ct.x_bin_center
ct = h2_metap_ct.y_bin_center

h2_metap_ct.xlabel = "$M(\pi^{+}\pi^{-}\eta)$"
h2_metap_ct.ylabel = "$cos\\theta_{GJ}$"
        
#%%

def fit_histo(h2d):
    
    A_a, x0_a, sigma_a, b0_a, c0_a  = [], [], [], [], []
    
    for i in range(h2d.nbins_y):
        '''
        h = h2d.project_x(bins = [i])
        M = h.bin_center
        C = h.bin_content
        dC = h.bin_error
    
        mr = np.linspace(M[0], M[-1], 1000)
        f = cf.gauss_fit()
        f.set_fit_list(fit = ['A', 'sigma',   'c0', 'b0'])
        f.fit_gaussbt(M, C, dC)
        A_a.append([f.A.value, f.A.err])
        sigma_a.append([f.sigma.value, f.sigma.err])
        b0_a.append([f.b0.value, f.b0.err])
        c0_a.append([f.c0.value, f.c0.err])
        plt.figure()
        h.plot_exp()
        f.plot_fit()
        B.plot_line(mr, f.gauss(mr))
        B.plot_line(mr, f.bkg(mr))
        '''
    
        h = h2d.project_x(bins=[i])
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
        f.c0.set(500)

        # set bounds for the fit parameters
        f.A_min.set(0.); f.A_max.set(1e5)
        f.x0_min.set(0.95); f.x0_max.set(0.96)
        f.sigma_min.set(0.008); f.sigma_max.set(0.016)
        f.c0_min.set(0.); f.c0_max.set(1e5)
        f.b0_min.set(0.00); f.b0_max.set(0.20)


        f.set_fit_list(fit = ['A', 'x0', 'sigma',  'c0', 'b0'])
        f.fit_gaussbt(M, C, dC) # gauss peak with bernstein bkg
        A_a.append([f.A.value, f.A.err])
        x0_a.append([f.x0.value, f.x0.err])
        sigma_a.append([f.sigma.value, f.sigma.err])
        b0_a.append([f.b0.value, f.b0.err])
        c0_a.append([f.c0.value, f.c0.err])
        plt.figure()
        h.plot_exp()
        f.plot_fit()
        B.plot_line(mr, f.gauss(mr))
        B.plot_line(mr, f.bt_bkg(mr))
        
        
    return    np.array(A_a), np.array(x0_a),  np.array(sigma_a), np.array(b0_a), np.array(c0_a)



#%%

A_a, x0_a, sigma_a, b0_a, c0_a  = fit_histo(h2_metap_ct)  

#%%

A_value = A_a[:,0] ; A_err = A_a[:,1] 
x0_value = x0_a[:,0] ; x0_err = x0_a[:,1] 
sigma_value = sigma_a[:,0] ; sigma_err = sigma_a[:,1] 
b0_value = b0_a[:,0]  ; b0_err = b0_a[:,1] 
c0_value = c0_a[:,0]  ; c0_err = c0_a[:,1] 

plt.figure();B.plot_exp(ct,  A_value, A_err,
                        x_label = "$cos\\theta_{GJ}$",  y_label = "A" )
plt.figure();B.plot_exp(ct,  x0_value, x0_err,
                        x_label = "$cos\\theta_{GJ}$",  y_label = "x0" )
plt.figure();B.plot_exp(ct,  sigma_value, sigma_err,
                        x_label = "$cos\\theta_{GJ}$",  y_label = "sigma" )
plt.figure();B.plot_exp(ct,  b0_value, b0_err, 
                        x_label = "$cos\\theta_{GJ}$",  y_label = "b0" )
plt.figure();B.plot_exp(ct,  c0_value, c0_err ,
                        x_label = "$cos\\theta_{GJ}$",  y_label = "c0" )
        
