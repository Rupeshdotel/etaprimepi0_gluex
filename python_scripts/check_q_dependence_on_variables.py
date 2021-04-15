#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:44:29 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import  matplotlib.pyplot  as plt
#%%

A = B.Parameter(500., 'A')
x0 = B.Parameter(0.956, 'x0')
sig = B.Parameter(.005, 'sigma')

def gaus(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2))


a0 = B.Parameter(1., 'a0')
a1 = B.Parameter(1., 'a1')


def lin_bkg(x):
    return a0() + a1()*x 

def signal(x):
    return gaus(x) + lin_bkg(x)



#%%
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18.npz')
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/variables_test_for_qfactor_be_t_phi.npz')

metap = d['metap']
metappi0 = d['metappi0']
cost_etap_gj = d['cost_etap_gj']
phi_etap_gj = d['phi_etap_gj']
mant = d['mant']
e_beam = d['e_beam']

#%%

cost_min = -1.0
cost_max = 1.0
bins_cost = 15

phi_min = -180.
phi_max = 180.
bins_phi = 15

t_min = 0.1
t_max = 0.7
bins_t = 15

be_min = 8.2
be_max = 8.8
bins_be = 15

metap_min = 0.85
metap_max = 1.05
metap_binwidth = 7e-3 # bin width in GeV
bins_metap = int((metap_max - metap_min)/metap_binwidth)





h2_etap_cost = B.histo2d(metap, cost_etap_gj,  bins = (bins_metap, bins_cost), 
                            range = [[metap_min, metap_max],[cost_min, cost_max]],
                       title = "Angular dependence of etaprime signal",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$\cos\\theta_{GJ}$" )

h2_etap_phi = B.histo2d(metap, phi_etap_gj,  bins = (bins_metap, bins_phi), 
                            range = [[metap_min, metap_max],[phi_min, phi_max]],
                       title = "Angular dependence of etaprime signal",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$\\phi_{GJ}$" )


h2_etap_mant = B.histo2d(metap, mant,  bins = (bins_metap, bins_t), 
                            range = [[metap_min, metap_max],[t_min, t_max]],
                       title = "t dependence of etaprime signal",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "momentun trasfer t" )

h2_etap_be = B.histo2d(metap, e_beam,  bins = (bins_metap, bins_be), 
                            range = [[metap_min, metap_max],[be_min, be_max]],
                       title = "beam energy dependence of etaprime signal",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "beam energy" )


metappi0_min = 1.0
metappi0_max = 2.5
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max -metappi0_min)/metappi0_binwidth)   

h2_ep_cost = B.histo2d(metappi0, cost_etap_gj, bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$\cos\\theta_{GJ}$" )

#%%

def fit_draw_histo(h2d, Fit = False):
    
    A_a, x0_a,  sigma_a, a0_a, a1_a, Q = [], [], [], [], [], []
    
    for i in range(h2d.nbins_y):
        h_xproj = h2d.project_x(bins=[i])
        M = h_xproj.bin_center
        
        
        minimum = 0.90
        maximum = 1.04
        sel  = (minimum < h_xproj.bin_center) & (h_xproj.bin_center < maximum)
        
        
        fit = B.genfit(signal, [A, x0, sig,  a0, a1],
                       x = h_xproj.bin_center[sel] , y = h_xproj.bin_content[sel], y_err = h_xproj.bin_error[sel],
                      print_results = False )
        
        
        #get the fit parameters
        #amplitude
        A_a.append([fit.parameters[0].value, fit.parameters[0].err])
        #mean
        x0_a.append([fit.parameters[1].value, fit.parameters[1].value])
        #sigma
        sigma_a.append([fit.parameters[2].value, fit.parameters[2].err])
        #intercept
        a0_a.append([fit.parameters[3].value, fit.parameters[3].err])
        #slope
        a1_a.append([fit.parameters[4].value, fit.parameters[4].err])
        
        #plot the data and the fit, both peak and bkg
        plt.figure()
        h_xproj.plot_exp()
        
        
        mr = np.linspace(M[sel].min(), M[sel].max(), 1000)
        B.plot_line(fit.xpl, fit.ypl)
        B.plot_line(mr, gaus(mr))
        B.plot_line(mr, lin_bkg(mr))
        
        ws = gaus(fit.parameters[1].value)
        wb = lin_bkg(fit.parameters[1].value)
        # calculate q-factor  
        q = ws/(ws+wb)
        Q.append(q)
            

        
    return  np.array(A_a), np.array(x0_a),  np.array(sigma_a), np.array(a0_a), np.array(a1_a), np.array(Q)
        
        
#%%

A_a,X_a,S_a,A0_a, A1_a, Q  = fit_draw_histo(h2_etap_phi, Fit = True)

#%%
A_a,X_a,S_a,A0_a, A1_a, Q  = fit_draw_histo(h2_etap_mant, Fit = True)


#%%
A_a,X_a,S_a,A0_a, A1_a, Q  = fit_draw_histo(h2_etap_be, Fit = True)
 
 
#%%
phi_bin_center = h2_etap_phi.y_bin_center

B.plot_exp(phi_bin_center, Q, marker = 'o')
B.pl.hlines(Q.mean(), phi_bin_center.min(), phi_bin_center.max())    
plt.xlabel('azimuthal angle $\phi_{GJ}$', fontsize=8)
plt.ylabel('Q values', fontsize=8)

#%%
mant_bin_center = h2_etap_mant.y_bin_center
B.plot_exp(mant_bin_center, Q, marker = 'o')
B.pl.hlines(Q.mean(), mant_bin_center.min(), mant_bin_center.max())    
plt.xlabel('momentum transfer t', fontsize=8)
plt.ylabel('Q values', fontsize=8)

#%%

e_beam_bin_center = h2_etap_be.y_bin_center
B.plot_exp(e_beam_bin_center, Q, marker = 'o')
B.pl.hlines(Q.mean(), e_beam_bin_center.min(),e_beam_bin_center.max())    
plt.xlabel('beam energy', fontsize=8)
plt.ylabel('Q values', fontsize=8)

#%%

mean_q = Q.mean()
se_q = Q.std()/np.sqrt(Q.shape)

error_percentage = se_q/mean_q *100

print(f' Mean of Qvalues = {mean_q:.4f}, standard error on qvalues = {se_q[0]:.4f} , error percentage = {error_percentage[0]:.4f}')



