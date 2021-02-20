#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:20:23 2020

@author: rupeshdotel
"""

#from mpl_toolkits import mplot3d
import numpy as np
import LT.box as B
import matplotlib.pyplot as plt
import class_fit as cf
#from class_fit import gauss_bt_fit

#%%
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_eventsprompt_gluexI.npz')


metap = d['etaprimemass_unique']
cos_GJ_etap = d['etaprimecosthetaGJ_unique']
mpi0 = d['pi0mass_unique']

pi0_min = 0.12
pi0_max = 0.15

etap_min = 0.86
etap_max = 1.03

sel = (pi0_min <  mpi0) &  (mpi0  < pi0_max) & (etap_min < metap) & (metap < etap_max)
xbin = 18
ybin = 18

h2 = B.histo2d( metap[sel], cos_GJ_etap[sel], bins = (xbin,ybin),
                   range=np.array([ (0.86, 1.03), (-1.0, 1.0)]), title = '2D')
epm = h2.x_bin_center
ct = h2.y_bin_center

h2.xlabel = "$M(\pi^{+}\pi^{-}\eta)$"
h2.ylabel = "$cos\\theta_{GJ}$"

#%%

def fit_histo(y_bin = 18):
    
    A_a, sigma_a, c0_a, b0_a = [], [], [], []
    
    for i in range(y_bin):
        
        h = h2.project_x(bins = [i])
        M = h.bin_center
        C = h.bin_content
        dC = h.bin_error
    
        mr = np.linspace(M[0], M[-1], 1000)
        f = cf.gauss_fit()
        f.set_fit_list(fit = ['A', 'sigma',   'c0', 'b0'])
        f.fit(M, C, dC)
        A_a.append([f.A.value, f.A.err])
        sigma_a.append([f.sigma.value, f.sigma.err])
        b0_a.append([f.b0.value, f.b0.err])
        c0_a.append([f.c0.value, f.c0.err])
        plt.figure()
        h.plot_exp()
        f.plot_fit()
        B.plot_line(mr, f.gauss(mr))
        B.plot_line(mr, f.bkg(mr))
        
    return    np.array(A_a),  np.array(sigma_a), np.array(b0_a), np.array(c0_a)
    


#%%

f = cf.gauss_fit()
x = np.linspace(0,1,1000)
#plt.plot(x, c.B23(x))
plt.plot(x, f.B34(x))

#%%




