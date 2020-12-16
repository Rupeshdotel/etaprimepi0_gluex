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

mep_bins = 18 # bins etaprime invariant mass
mp_bins = 18 # bins pi0 invariant mass
mep_min = 0.86 # left edge of etaprime invariant mass
mep_max = 1.06 # right edge of etaprime invariant mass
mp_min = 0.12 # left edge of pi0 invariant mass
mp_max = 0.15 # right edge of pi0 invariant mass
h_epi0 = B.histo2d( metap, mpi0, bins = (mep_bins, mp_bins),
                   range=np.array([ (mep_min, mep_max), (mp_min, mp_max)]), title = '2D_prompt', xlabel = "$M(\pi^{+}\pi^{-}\eta)$",
                   ylabel = "$M(\gamma\gamma)$")


epm = h_epi0.x_bin_center
pi0m = h_epi0.y_bin_center

def fit_histo(h2d):
     #define empty lists for parameters and their errors from  fits
     
     A_a, x0_a,  sigma_a, c0_a, b0_a = [], [], [], [], []
     

     for i in range(h2d.nbins_y): # looping over pi0mass bins
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
         f.b0.set(0.14)
         f.c0.set(50000)
         
         # set bounds for the fit parameters
         f.A_min.set(0.); f.A_max.set(1e5)
         f.x0_min.set(0.94); f.x0_max.set(0.98)
         f.sigma_min.set(0.008); f.sigma_max.set(0.016)
         f.c0_min.set(0.); f.c0_max.set(1e5)
         f.b0_min.set(0.10); f.b0_max.set(0.20)
         
         
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
A_a, x0_a, S_a, B0_a, C0_a = fit_histo(h_epi0)    

#%%
A_value = A_a[:,0] 
A_err = A_a[:,1] 

x0_value = x0_a[:,0] 
x0_err = x0_a[:,1] 

S_value = S_a[:,0] 
S_err = S_a[:,1] 

B0_value = B0_a[:,0] 
B0_err = B0_a[:,1] 

C0_value = C0_a[:,0] 
C0_err = C0_a[:,1] 


e = cf.gauss_fit()
e.set_fit_list( fit = ['A', 'x0', 'sigma', 'k0', 'k1'])

e.A.set(2500)
e.x0.set(0.135)
e.sigma.set(0.005)
e.k0.set(0.)
e.k1.set(0.)

e.A_min.set(1000); e.A_max.set(3000)
e.x0_min.set(0.13); e.x0_max.set(0.14)
e.sigma_min.set(0.0001); e.sigma_max.set(0.01)

e.k0_min.set(-1e5); e.k0_max.set(1e5)
e.k1_min.set(-1e5); e.k1_max.set(1e5)

e.fit_gausslin(pi0m, A_value, A_err) # gauss peak with linear bkg
figure();e.plot_fit()
B.plot_line(pi0m, e.gauss(pi0m))
B.plot_line(pi0m, e.lin_bkg(pi0m))
B.plot_exp(pi0m, A_a[:,0], A_a[:,1],  plot_title = 'Fit the fit parameter A',  x_label = ' $M(\gamma\gamma)$' )


plt.figure()
px0 = B.polyfit(pi0m, x0_value, x0_err,   order = 1)
B.plot_exp(pi0m, x0_value,  x0_err, plot_title = 'Fit the fit parameter $x0$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(px0.xpl, px0.ypl)

plt.figure()
pS = B.polyfit(pi0m, S_value, S_err,   order = 2)
B.plot_exp(pi0m, S_value,  S_err, plot_title = 'Fit the fit parameter $sigma$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pS.xpl, pS.ypl)


plt.figure()
pB0 = B.polyfit(pi0m, B0_value, B0_err,   order = 2)
B.plot_exp(pi0m, B0_value,  B0_err, plot_title = 'Fit the fit parameter $b0$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pB0.xpl, pB0.ypl)



plt.figure()
pC0 = B.polyfit(pi0m, C0_value, C0_err,   order = 4)
B.plot_exp(pi0m, C0_value,  C0_err, plot_title = 'Fit the fit parameter $c0$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
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
    mean = px0.poly(y)
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



sel = (mp_min <  mpi0) &  (mpi0  < mp_max) & (mep_min < metap) & (metap < mep_max)

#x = etaprimepi0mass_unique_sc[sel]
#y = etaprimecosthetaGJ_unique_sc[sel]
qb_sel = qb[sel]
qs_sel = qs[sel]

#%%

hb_ep = B.histo(metap[sel], range = (mep_min, mep_max),  bins = 24, weights = qb_sel, title = "bkg $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
hs_ep = B.histo(metap[sel], range = (mep_min, mep_max), bins = 24, weights = qs_sel, title = "signal $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
ht_ep = B.histo(metap[sel], range = (mep_min, mep_max), bins = 24, title = 'Non-weighted',
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

# for pi0mass
hb_p = B.histo(mpi0[sel], range = (mp_min, mp_max),  bins = 24, weights = qb_sel, title = 'bkg $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
hs_p = B.histo(mpi0[sel], range = (mp_min, mp_max), bins = 24, weights = qs_sel, title = 'signal $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
ht_p = B.histo(mpi0[sel], range = (mp_min, mp_max), bins = 24, title = 'Non-weighted',
               xlabel = "$M(\gamma\gamma)$")
    
hs2 = B.histo2d(metap[sel], mpi0[sel], range = [[mep_min, mep_max], [mp_min, mp_max]],
                bins = (24, 24), title = 'Invariant Mass in 2D for Signal events', weights = qs_sel, 
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
hb2 = B.histo2d(metap[sel], mpi0[sel], range = [[mep_min, mep_max], [mp_min, mp_max]],
                bins = (24, 24), title = 'Invariant Mass in 2D for Bkg events', weights = qb_sel,
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
ht2 = B.histo2d(metap[sel], mpi0[sel], range = [[mep_min, mep_max], [mp_min, mp_max]], bins = (24, 24), 
                title = 'Invariant Mass in 2D ', xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
    
    
#%%
# for checking 
i = 4
figure();
h_epi0.project_x(bins = [i]).plot_exp()
F = cf.gauss_fit()
F.A.set(A_a[i][0])
F.x0.set(x0_a[i][0])
F.sigma.set(S_a[i][0])
F.b0.set(B0_a[i][0])
F.c0.set(C0_a[i][0])
mr = np.linspace(mep_min, mep_max, 1000)
plt.plot(mr, F.signal_bt_bkg(mr))





