#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 18:56:30 2020

@author: rupeshdotel
"""

import numpy as np
import matplotlib.pyplot as plt
import LT.box as B

#%%
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_eventsprompt_gluexI.npz')

mpi0 = d['pi0mass_unique']
metap = d['etaprimemass_unique']



phi_GJ_pi0 = d['pi0costhetaGJ_unique']
phi_GJ_etap = d['etaprimecosthetaGJ_unique']


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
q = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qvalues_gluexI.npz')

qs_sel = q['qs_sel']; qb_sel = q['qb_sel'];
sel=q['sel']; 

#%%

he_ws2d = B.histo(meta[sel], bins = 40, weights = qs_sel, title = '$\eta$' , xlabel = 'M($\gamma\gamma$)')
he = B.histo(meta[sel], bins = 40)
he_ws2d.plot_exp()
he_ws2d.show_fit_list()
he_ws2d.set_fit_list(['A', 'mean', 'sigma', 'b0', 'b1'])
he_ws2d.fit()
he_ws2d.plot_fit()
he_ws2d.title = '$\eta$'
he_ws2d.x_lable = '$M(\gamma\gamma)$'
b0 = he_ws2d.b0.value
b1 = he_ws2d.b1.value
#b2 = hs.b2.value

def eta_bkg(x):
    return b0 + b1*x 

def qeta_bkg(x):
    b = eta_bkg(x)
    t = he_ws2d.fit_func(x) 
    s = t - b
    qs = s/t
    qb = 1 - qs
    return qs, qb



qse, qbe = qeta_bkg(meta[sel])

hse = B.histo(meta[sel], bins = 40, weights = qs_sel*qse, title = '$\eta$', xlabel = 'M($\gamma\gamma$)')
plt.figure();he.plot_exp(); hse.plot_exp();  he_ws2d.plot_exp() #

#%%

hep = B.histo(metap[sel], bins = 30, title = "$\eta^{'}$" , xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

A = B.Parameter(17500., 'A') #gluexI
#A = B.Parameter(2500., 'A')
x0 = B.Parameter(0.955, 'x0')
s = B.Parameter(.01, 's')



def gaus_peak(x):
    return A() * np.exp(-(x - x0())**2/(2.*s()**2))


alpha = B.Parameter(3000., 'alpha')
beta = B.Parameter(3000., 'beta')
gamma = B.Parameter(10., 'gamma')
delta = B.Parameter(1., 'delta')

def tanh_bkg(x):
    return alpha() + beta()*np.tanh(gamma()*(x - delta()))

def signal(x):
    return tanh_bkg(x) + gaus_peak(x)


M = hep.bin_center
C = hep.bin_content
dC = hep.bin_error
fit = B.genfit(signal, [ A, x0, s, alpha, beta, gamma, delta], x = M, y = C, y_err = dC )

plt.figure()
B.plot_exp(M, C, dC, x_label = '$M(\pi^{+}\pi^{-}\eta)$', plot_title = 'Gaussian peak on Tanh bkg   ')

B.plot_line(fit.xpl, fit.ypl)
             


fit.save_parameters()
al = fit.parameters_sav[3].value
be = fit.parameters_sav[4].value
ga = fit.parameters_sav[5].value
de = fit.parameters_sav[6].value


def etap_bkg(x):
    return al + be*np.tanh(ga*(x - de))


 

def qetap_bkg(x):
    t = fit.func(x)
    b = etap_bkg(x)
    s = t - b
    qs = s/t
    qb = 1 - qs
    return qs, qb



qsep, qbep = qetap_bkg(metap[sel])
qsep = np.abs(qsep)
qbep = np.abs(qbep)

hsep = B.histo(metap[sel], bins = 30, weights = qsep, title = "$\eta^{'}$" , xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
hbep = B.histo(metap[sel], bins = 30, weights = qbep, title = "$\eta^{'}$" , xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
hep = B.histo(metap[sel], bins = 30)

plt.figure();
hep.plot_exp();  hbep.plot_exp(); hsep.plot_exp();


#%%

#hp = B.histo(mpi0[sel], bins = 30, weights = qs_sel)
hp = B.histo(mpi0[sel], bins = 18,  title = '$\pi^{0}$' , xlabel = 'M($\gamma\gamma$)')
plt.figure()
hp.plot_exp()
hp.show_fit_list()
hp.set_fit_list(['A', 'mean', 'sigma', 'b0', 'b1'])
hp.fit()
hp.plot_fit()
b0 = hp.b0.value
b1 = hp.b1.value
#b2 = hs.b2.value

def pi0_bkg(x):
    return b0 + b1*x 

def qpi0_bkg(x):
    b = pi0_bkg(x)
    t = hp.fit_func(x) 
    s = t - b
    qs = s/t
    qb = 1 - qs
    return qs, qb



qsp, qbp = qpi0_bkg(mpi0[sel])
qsp = np.abs(qsp)
qbp = np.abs(qbp)



hsp = B.histo(mpi0[sel], bins = 18, weights = qsp, title = '$\pi^{0}$' , xlabel = 'M($\gamma\gamma$)')
hsb = B.histo(mpi0[sel], bins = 18, weights = qbp, title = '$\pi^{0}$' , xlabel = 'M($\gamma\gamma$)')
plt.figure(); hp.plot_exp(); hsb.plot_exp(); hsp.plot_exp();  


#%%
#for etaprime mass
hb_ep = B.histo(metap[sel], range = (0.86, 1.06),  bins = 24, weights = qbep, title = "bkg $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
hs_ep = B.histo(metap[sel], range = (0.86, 1.06), bins = 24, weights = qsep*qse, title =  "signal $\eta'$",
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")
ht_ep = B.histo(metap[sel], range = (0.86, 1.06), bins = 24, title = 'Non-weighted',
               xlabel = "$M(\pi^{+}\pi^{-}\eta)$")

# for pi0mass
hb_p = B.histo(mpi0[sel], range = (0.12, 0.15),  bins = 12, weights = qbp, title = 'bkg $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
hs_p = B.histo(mpi0[sel], range = (0.12, 0.15), bins = 12, weights = qsp, title = 'signal $\pi^{0}$',
               xlabel = "$M(\gamma\gamma)$")
ht_p = B.histo(mpi0[sel], range = (0.12, 0.15), bins = 12, title = 'Non-weighted',
               xlabel = "$M(\gamma\gamma)$")


# for etamass
hb_eta = B.histo(meta[sel], range = (0.48, 0.60),  bins = 50, weights = qbe, title = 'bkg',
               xlabel = "$M(\gamma\gamma)$")
hs_eta = B.histo(meta[sel], range = (0.48, 0.60), bins = 50, weights = qse, title = 'sig',
               xlabel = "$M(\gamma\gamma)$")
ht_eta = B.histo(meta[sel], range = (0.48, 0.60), bins = 50, title = 'All',
               xlabel = "$M(\gamma\gamma)$")




hs2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.06], [0.12, 0.15]],
                bins = (20, 10), title = '2d_sig', weights = qsep*qsp*qse, 
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
hb2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.06], [0.12, 0.15]],
                bins = (20, 10), title = '2d_bkg', weights = qbep*qbp,
                xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")
ht2 = B.histo2d(metap[sel], mpi0[sel], range = [[0.86, 1.06], [0.12, 0.15]], bins = (20, 10), 
                title = '2d_total', xlabel = "$M(\pi^{+}\pi^{-}\eta)$", ylabel = "$M(\gamma\gamma)$")


hs2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.3, 3], [-1, 1]],
                bins = (50, 40), title = '2d_sig_X_GJ', weights = qsep*qsp*qse, 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

hb2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.3, 3], [-1, 1]],
                bins = (50, 40), title = '2d_bkg_X_GJ', weights = qbep*qbp, 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

ht2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.3, 3], [-1, 1]],
                bins = (50, 40), title = '2d_X_GJ',
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")








#%%
hs2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.3, 3], [-1, 1]],
                bins = (50, 40), title = '2d_sig_X_GJ', weights = qsep*qsp*qse, 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

hb2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.3, 3], [-1, 1]],
                bins = (50, 40), title = '2d_bkg_X_GJ', weights = qbep*qbp, 
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")

ht2_X_GJ = B.histo2d(metap_pi0[sel], cos_GJ_etap[sel], range = [[1.3, 3], [-1, 1]],
                bins = (50, 40), title = '2d_X_GJ',
                xlabel = "$M(\eta^{'}\pi^{0})$", ylabel = "$cos\\theta_{GJ}$")



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













