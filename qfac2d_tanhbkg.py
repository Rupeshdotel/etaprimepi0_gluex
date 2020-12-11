#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 12:40:43 2020

@author: rupeshdotel
"""
from mpl_toolkits import mplot3d
import numpy as np
import LT.box as B
import matplotlib.pyplot as plt

#%%

Ae = B.Parameter(2500., 'Ae')
x0e = B.Parameter(0.956, 'x0e')
se = B.Parameter(.01, 'se')

alpha = B.Parameter(10000., 'alpha')
beta = B.Parameter(10000., 'beta')
gamma = B.Parameter(10., 'gamma')
delta = B.Parameter(1., 'delta')


'''#S17
Ae = B.Parameter(150., 'Ae')
x0e = B.Parameter(0.956, 'x0e')
se = B.Parameter(.01, 'se')
alpha = B.Parameter(50., 'alpha')
beta = B.Parameter(100., 'beta')
gamma = B.Parameter(10., 'gamma')
delta = B.Parameter(1., 'delta')'''


def gaus_peak(x):
    return Ae() * np.exp(-(x - x0e())**2/(2.*se()**2))



def tanh_bkg(x):
    return alpha() + beta()*np.tanh(gamma()*(x - delta()))

def signal(x):
    return tanh_bkg(x) + gaus_peak(x)

#B.plot_exp(x, signal(x))

#%%

#file = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_eventsprompt_17.npz')

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
xbin = 30
ybin = 18
h_epi0 = B.histo2d( metap, mpi0, bins = (xbin,ybin),
                   range=np.array([ (0.86, 1.06), (0.12, 0.15)]), title = '2D_prompt', xlabel = "$M(\pi^{+}\pi^{-}\eta)$",
                   ylabel = "$M(\gamma\gamma)$")


epm = h_epi0.x_bin_center
pi0m = h_epi0.y_bin_center

def fit_histo(epm = True):
     #define empty lists for parameters and their errors from  fits
     
     Aa = []; Aa_err = []
     X0 = []; X0_err = []
     S = []; S_err = []
     
     Al = []; Al_err = []
     Be = []; Be_err = []
     Ga = []; Ga_err = []
     De = []; De_err = []
     

     for i in range(h_epi0.bin_content.shape[1]):
         h = h_epi0.project_x(bins=[i])
         M = h.bin_center
         C = h.bin_content
         dC = h.bin_error
         fit = B.genfit(signal, [ Ae, x0e, se, alpha, beta, gamma, delta], x = M, y = C, y_err = dC )
         
         if i == 9:
             plt.figure()
             #B.plot_exp(M, C, dC, x_label = '$M(\pi^{+}\pi^{-}\eta)$', plot_title = 'Gaussian peak on Tanh bkg  ')
             h.plot_exp()
             B.plot_line(fit.xpl, fit.ypl)
             
             
             #fit.parameters_sav[0].value
             #fit.parameters_sav[0].err
         fit.save_parameters()
         Aa.append(fit.fit_result.x[0]); Aa_err.append(fit.parameters_sav[0].err)
         X0.append(fit.fit_result.x[1]); X0_err.append(fit.parameters_sav[1].err)
         S.append(fit.fit_result.x[2]);  S_err.append(fit.parameters_sav[2].err)
         Al.append(fit.fit_result.x[3]); Al_err.append(fit.parameters_sav[3].err)
         Be.append(fit.fit_result.x[4]); Be_err.append(fit.parameters_sav[4].err)
         Ga.append(fit.fit_result.x[5]); Ga_err.append(fit.parameters_sav[5].err)
         De.append(fit.fit_result.x[6]); De_err.append(fit.parameters_sav[6].err)
         
     Aa = np.array(Aa); Aa_err = np.array(Aa_err)
     X0 = np.array(X0); X0_err = np.array(X0_err)
     S = np.array(S); S_err = np.array(S_err)
     Al = np.array(Al); Al_err = np.array(Al_err)
     Be = np.array(Be); Be_err = np.array(Be_err)
     Ga = np.array(Ga); Ga_err = np.array(Ga_err)
     De = np.array(De); De_err = np.array(De_err)
         
    
     
        
     return Aa, X0, S, Al, Be, Ga, De, Aa_err, X0_err, S_err, Al_err, Be_err, Ga_err, De_err
#%% 
Aa, X0, S, Al, Be, Ga, De,  Aa_err, X0_err, S_err, Al_err, Be_err, Ga_err, De_err = fit_histo()

#%%
'''pAl = B.polyfit(pi0m, Al, Al_err,  order=4)
pBe = B.polyfit(pi0m, Be, Be_err, order=4)
pGa = B.polyfit(pi0m, Ga, Ga_err, order=4)
pDe = B.polyfit(pi0m, De, De_err,  order=3)


pA = B.polyfit(pi0m, Aa, Aa_err,  order=4)
pS = B.polyfit(pi0m, S, S_err, order=1)
pM = B.polyfit(pi0m, X0, X0_err, order=1)

pAl = B.polyfit(pi0m, Al, Al_err,  order=4)
B.plot_exp(pi0m, Al, Al_err, plot_title = "Fit the fit parameter $\\alpha$", x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pAl.xpl, pAl.ypl)

plt.figure()
pBe = B.polyfit(pi0m, Be, Be_err,  order=4)
B.plot_exp(pi0m, Be, Be_err, plot_title = "Fit the fit parameter $\\beta$", x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pBe.xpl, pBe.ypl)


plt.figure()
pAa = B.polyfit(pi0m, Aa, Aa_err,  order=4)
B.plot_exp(pi0m, Aa, Aa_err, plot_title = "Fit the fit parameter $height A$", x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pAa.xpl, pAa.ypl)
'''

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


fit_Aa = B.genfit(signal_p, [A, x0, sig, a0, a1], x = pi0m, y = Aa, y_err = Aa_err )
plt.figure()
B.plot_line(fit_Aa.xpl, fit_Aa.ypl)
#B.plot_line(pi0m, gaus_p(pi0m))
#B.plot_line(pi0m, lin_bkgp(pi0m))
B.plot_exp(pi0m, Aa, Aa_err,  plot_title = 'Fit the fit parameter A',  x_label = ' $M(\gamma\gamma)$' )


plt.figure()
pM = B.polyfit(pi0m, X0, X0_err,  order=1)
B.plot_exp(pi0m, X0, X0_err, plot_title = 'Fit the fit parameter $mean$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pM.xpl, pM.ypl)

plt.figure()
pS = B.polyfit(pi0m, S, S_err,  order=1)
B.plot_exp(pi0m, S, S_err, plot_title = 'Fit the fit parameter $sigma$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pS.xpl, pS.ypl)

fit_al = B.genfit(signal_p, [A, x0, sig, a0, a1], x = pi0m, y = Al, y_err = Al_err )
plt.figure()
B.plot_line(fit_al.xpl, fit_al.ypl)
B.plot_line(pi0m, gaus_p(pi0m))
B.plot_line(pi0m, lin_bkgp(pi0m))
B.plot_exp(pi0m, Al, Al_err,  plot_title = 'Fit the fit parameter $\\alpha$',  x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$' )
#B.plot_exp(pi0m, Al, Al_err)


fit_be = B.genfit(signal_p, [A, x0, sig, a0, a1], x = pi0m, y = Be, y_err = Be_err )
plt.figure()
B.plot_line(fit_be.xpl, fit_be.ypl)
B.plot_line(pi0m, gaus_p(pi0m))
B.plot_line(pi0m, lin_bkgp(pi0m))
B.plot_exp(pi0m, Be, Be_err,  plot_title = 'Fit the fit parameter $\\beta$',  x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$' )
#B.plot_exp(pi0m, Be, Be_err)

plt.figure()
pGa = B.polyfit(pi0m, Ga, Ga_err,  order=4)
B.plot_exp(pi0m, Ga, Ga_err, plot_title = 'Fit the fit parameter $\gamma$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pGa.xpl, pGa.ypl)

plt.figure()
pDe = B.polyfit(pi0m, De, De_err,  order=1)

B.plot_exp(pi0m, De, De_err, plot_title = 'Fit the fit parameter $\delta$', x_label = 'bin centers of y-axis in 2D plot $M(\pi^{0})$')
B.plot_line(pDe.xpl, pDe.ypl)





#%%
#construct 
'''def bkg_fit(x,y):
    y = pAl.poly(y) + pBe.poly(y) * np.tanh(pGa.poly(y) * (x - pDe.poly(y))) 
    return y

def peak_fit(x,y):
    A = pA.poly(y)
    mean = pM.poly(y)
    sigma = pS.poly(y)
    y = A*np.exp(-0.5*((x-mean)/sigma)**2)
    return y


def signal_fit(x,y):
    return peak_fit(x,y) + bkg_fit(x,y)'''

def bkg_fit(x,y):
    y = fit_al.func(y) + fit_be.func(y) * np.tanh(pGa.poly(y) * (x - pDe.poly(y))) 
    return y

def peak_fit(x,y):
    A = fit_Aa.func(y)
    mean = pM.poly(y)
    sigma = pS.poly(y)
    y = A*np.exp(-0.5*((x-mean)/sigma)**2)
    return y


def signal_fit(x,y):
    return peak_fit(x,y) + bkg_fit(x,y)

#%%
q_b = np.abs(bkg_fit(metap, mpi0))
q_s = np.abs(peak_fit(metap, mpi0))

q = q_s/(q_b + q_s)
qs = np.abs(q)
qb = 1-q

#%%
pi0_min = 0.12
pi0_max = 0.15

etap_min = 0.86
etap_max = 1.06

sel = (pi0_min <  mpi0) &  (mpi0  < pi0_max) & (etap_min < metap) & (metap < etap_max)

#x = etaprimepi0mass_unique_sc[sel]
#y = etaprimecosthetaGJ_unique_sc[sel]
qb_sel = qb[sel]
qs_sel = qs[sel]

qeta = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qvalues_gluexI.npz')
qse = qeta['qse']

           

#%%

#1d histos

#for etaprime mass
hb_ep = B.histo(metap[sel], range = (0.86, 1.06),  bins = xbin, weights = qb_sel, title = "bkg $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
hs_ep = B.histo(metap[sel], range = (0.86, 1.06), bins = xbin, weights = qs_sel*qse, title = "signal $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
ht_ep = B.histo(metap[sel], range = (0.86, 1.06), bins = xbin, title = 'Non-weighted',
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

# for pi0mass
hb_p = B.histo(mpi0[sel], range = (0.12, 0.15),  bins = ybin, weights = qb_sel, title = 'bkg $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
hs_p = B.histo(mpi0[sel], range = (0.12, 0.15), bins = ybin, weights = qs_sel*qse, title = 'signal $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
ht_p = B.histo(mpi0[sel], range = (0.12, 0.15), bins = ybin, title = 'Non-weighted',
               xlabel = "$M(\gamma\gamma)$")




hs2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.06], [0.12, 0.15]],
                bins = (xbin, ybin), title = 'Invariant Mass in 2D for Signal events', weights = qs_sel, 
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
hb2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.06], [0.12, 0.15]],
                bins = (xbin, ybin), title = 'Invariant Mass in 2D for Bkg events', weights = qb_sel,
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
ht2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.06], [0.12, 0.15]], bins = (xbin, ybin), 
                title = 'Invariant Mass in 2D ', xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")


hs2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.0, 3], [-1, 1]],
                bins = (50, 40), title = 'Angular Distribution in Gottfried-Jackson Frame for Signal Events', weights = qs_sel*qse, 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

hb2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.0, 3], [-1, 1]],
                bins = (50, 40), title = 'Angular Distribution in Gottfried-Jackson Frame for Bkg Events', weights = qb_sel, 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

ht2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.0, 3], [-1, 1]],
                bins = (50, 40), title = 'Angular Distribution in Gottfried-Jackson Frame ',
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

#%%
'''#%%
h = B.histo(metap, bins=50, range=(0.85, 1.06))
M = h.bin_center
C = h.bin_content
dC = h.bin_error

B.plot_exp(M, C, dC)

#fit = B.genfit(signal, [ A, x0, s, xm, a1, a2, a3, a4], x = M, y = C, y_err = dC )
fit = B.genfit(signal, [ A, x0, s, alpha, beta, gamma, delta], x = M, y = C, y_err = dC )

B.plot_exp(M, C, dC)
B.plot_line(fit.xpl, fit.ypl)

x = np.linspace(0.85,1.06,1000)

a1 = B.Parameter(1., 'a1')
a2 = B.Parameter(5., 'a2')
a3 = B.Parameter(50., 'a3')
a4 = B.Parameter(35000., 'a4')
xm = B.Parameter(.85, 'xm')

def pol_bkg(x):
    return ((x - xm()) *  (a1()*x + a2()*x**2 + a3()*x**3 + a4()*x**4 ))


'''
#%%

#cuts selection
omega_thr = 0.85
#cut_omega = ( mpippimpi0 > omega_thr)
omega_min = 0.75
omega_max = 0.85
cut_omega = (omega_min > mpippimpi0[sel]) |  (mpippimpi0[sel] > omega_max)



#baryon_thr = 1.35
#cut_deltap = mpi0p > baryon_thr
#cut_deltapp = mpipp > baryon_thr

#cuts =  cut_omega & cut_deltap & cut_deltapp
cuts = cut_omega

mpi0_rw = mpi0[sel][cuts]
metap_rw = metap[sel][cuts]


phi_GJ_pi0_rw = phi_GJ_pi0[sel][cuts]
phi_GJ_etap_rw = phi_GJ_etap[sel][cuts]


cos_GJ_pi0_rw = cos_GJ_pi0[sel][cuts]
cos_GJ_etap_rw = cos_GJ_etap[sel][cuts]


meta_rw = meta[sel][cuts]
metap_pi0_rw = metap_pi0[sel][cuts]
BE_rw = BE[sel][cuts]
mpippimpi0_rw = mpippimpi0[sel][cuts]
mpi0p_rw = mpi0p[sel][cuts]
mpipp_rw = mpipp[sel][cuts]
num_combos_rw = num_combos[sel][cuts]
combo_num_rw = combo_num[sel][cuts]
dt_rw = dt[sel][cuts]
event_num_rw = event_num[sel][cuts]
kinfit_CL_rw = kinfit_CL[sel][cuts]
chisq_ndf_rw = chisq_ndf[sel][cuts]

qs_sel_rw = qs_sel[cuts]
qse_rw = qse[cuts]
















