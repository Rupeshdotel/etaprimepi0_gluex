#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:26:54 2020

@author: rupeshdotel
"""

from mpl_toolkits import mplot3d
import numpy as np
import LT.box as B
import matplotlib.pyplot as plt
import class_fit as cf
#%%
#load the data
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_eventsprompt_gluexI.npz')
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_events_accid.npz')
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_eventsprompt_17.npz')
mpi0 = d['pi0mass_unique']
metap = d['etaprimemass_unique']



phi_GJ_pi0 = d['pi0phiGJ_unique']
phi_GJ_etap = d['etaprimephiGJ_unique']


cos_GJ_pi0 = d['pi0costhetaGJ_unique']
cos_GJ_etap = d['etaprimecosthetaGJ_unique']


meta = d['etamass_unique']
metap_pi0 = d['etaprimepi0mass_unique']
BE = d['BE_unique']
mpippimpi0 = d['pippimpi0_unique']
mpi0p = d['pi0p_unique']
mpipp = d['pipp_unique']
num_combos = d['num_combos_unique']
combo_num = d['combo_number_unique']
dt = d['dt_unique']
event_num = d['event_num_unique']
kinfit_CL = d['kinfit_CL_unique']
chisq_ndf = d['chisq_ndf_unique']
#%%

xbin = 18
ybin = 18
h_epi0 = B.histo2d( metap, mpi0, bins = (xbin,ybin),
                   range=np.array([ (0.86, 1.04), (0.12, 0.15)]), title = '2D_prompt', xlabel = "$M(\pi^{+}\pi^{-}\eta)$",
                   ylabel = "$M(\gamma\gamma)$")


epm = h_epi0.x_bin_center
pi0m = h_epi0.y_bin_center

def fit_histo(y_bin = 18):
     #define empty lists for parameters and their errors from  fits
     
     A_a, sigma_a, c0_a, b0_a = [], [], [], []
     

     for i in range(y_bin):
         h = h_epi0.project_x(bins=[i])
         M = h.bin_center
         C = h.bin_content
         dC = h.bin_error
         mr = np.linspace(M[0], M[-1], 1000)
         f = cf.gauss_bt_fit()
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
A_a, S_a, B0_a, C0_a = fit_histo()    

#%%
#h_al = B.plot_exp(pi0m, Al, Al_err)

A = B.Parameter(2500., 'A') #Spring2018 Fall2018
#A = B.Parameter(100., 'A') #Spring2017
x0 = B.Parameter(0.135, 'x0')
sig = B.Parameter(.005, 'sigma')

def gaus_p(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2))

a0 = B.Parameter(0., 'a0')
a1 = B.Parameter(0., 'a1')
#a2 = B.Parameter(0., 'a2')

def lin_bkgp(x):
    return a0() + a1()*x 


def signal_p(x):
    return gaus_p(x) + lin_bkgp(x)
    

A_value = A_a[:,0] 
A_err = A_a[:,1] 

S_value = S_a[:,0] 
S_err = S_a[:,1] 

B0_value = B0_a[:,0] 
B0_err = B0_a[:,1] 

C0_value = C0_a[:,0] 
C0_err = C0_a[:,1] 

fit_Aa = B.genfit(signal_p, [A, x0, sig, a0, a1], x = pi0m, y = A_value, y_err = A_err )
plt.figure()
B.plot_line(fit_Aa.xpl, fit_Aa.ypl)
B.plot_line(pi0m, gaus_p(pi0m))
B.plot_line(pi0m, lin_bkgp(pi0m))
B.plot_exp(pi0m, A_a[:,0], A_a[:,1],  plot_title = 'Fit the fit parameter A',  x_label = ' $M(\gamma\gamma)$' )


plt.figure()
pS = B.polyfit(pi0m, S_a[:,0], S_a[:,1],   order = 1)
B.plot_exp(pi0m, S_a[:,0],  S_a[:,1], plot_title = 'Fit the fit parameter $sigma$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pS.xpl, pS.ypl)


plt.figure()
pB0 = B.polyfit(pi0m, B0_a[:,0], B0_a[:,1],   order = 2)
B.plot_exp(pi0m, B0_a[:,0],  B0_a[:,1], plot_title = 'Fit the fit parameter $b0$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pB0.xpl, pB0.ypl)



plt.figure()
pC0 = B.polyfit(pi0m, C0_a[:,0], C0_a[:,1],   order = 4)
B.plot_exp(pi0m, C0_a[:,0],  C0_a[:,1], plot_title = 'Fit the fit parameter $c0$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pC0.xpl, pC0.ypl)

#%%
def bkg_fit(x,y):
    #f = cf.gauss_bt_fit()
    #y = f.b34(pC0.poly(y), x, 0.86, 1.04)
    b0 = pB0.poly(y)
    
    db0 = 0.50
    m1 = 1.04
    m0 = 0.86
    b1 = b0 + db0
    
    
    a = db0/(m1 - m0)
    #b = (self.b1 * self.m0() - self.b0() * self.m1())/(self.m0() - self.m1())
    b = (b1 * m0 - b0 * m1) / (m0 - m1)
    X = a * x + b
    Y = 4 * X**3 * (1 - X)
    return   pC0.poly(y)  * Y
    

def peak_fit(x,y):
    A = fit_Aa.func(y)
    mean = 0.956
    sigma = pS.poly(y)
    y = A*np.exp(-0.5*((x-mean)/sigma)**2)
    return y


def signal_fit(x,y):
    return peak_fit(x,y) + bkg_fit(x,y)

#%%
    
q_b = np.abs(bkg_fit(metap, mpi0))
q_s = np.abs(peak_fit(metap, mpi0))


q = q_s/(q_b + q_s)
qs = q
qb = 1-q

pi0_min = 0.12
pi0_max = 0.15

etap_min = 0.86
etap_max = 1.04

sel = (pi0_min <  mpi0) &  (mpi0  < pi0_max) & (etap_min < metap) & (metap < etap_max)

#x = etaprimepi0mass_unique_sc[sel]
#y = etaprimecosthetaGJ_unique_sc[sel]
qb_sel = qb[sel]
qs_sel = qs[sel]

#%%

hb_ep = B.histo(metap[sel], range = (0.86, 1.04),  bins = 24, weights = qb_sel, title = "bkg $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
hs_ep = B.histo(metap[sel], range = (0.86, 1.04), bins = 24, weights = qs_sel, title = "signal $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
ht_ep = B.histo(metap[sel], range = (0.86, 1.04), bins = 24, title = 'Non-weighted',
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

# for pi0mass
hb_p = B.histo(mpi0[sel], range = (0.12, 0.15),  bins = 24, weights = qb_sel, title = 'bkg $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
hs_p = B.histo(mpi0[sel], range = (0.12, 0.15), bins = ybin, weights = qs_sel, title = 'signal $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
ht_p = B.histo(mpi0[sel], range = (0.12, 0.15), bins = ybin, title = 'Non-weighted',
               xlabel = "$M(\gamma\gamma)$")
    
hs2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.04], [0.12, 0.15]],
                bins = (24, 24), title = 'Invariant Mass in 2D for Signal events', weights = qs_sel, 
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
hb2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.04], [0.12, 0.15]],
                bins = (24, 24), title = 'Invariant Mass in 2D for Bkg events', weights = qb_sel,
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
ht2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.04], [0.12, 0.15]], bins = (24, 24), 
                title = 'Invariant Mass in 2D ', xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
    
    
