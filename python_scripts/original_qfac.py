#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 08:29:09 2021

@author: rupeshdotel
"""



import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import LT.box as B



#%%


#normalized gaussian
def p_peak(x,  mu,  sig, A):
     p = 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-0.5*((x-mu)/(1.*sig)**2))
     return A * p
     

    


# lin. background probability dist. (normalized)
def p_bkg(x, x0, x1, y0, y1):
    N = 0.5*(y0+y1)*(x1-x0)
    p = 1./N*((y1-y0)/(x1-x0)*(x-x0) + y0)
    return p


# combined signal probability dist.
def p_signal(x, P, x0, x1):
    fs = P[0]
    y0 = P[1]
    y1 = P[2]
    mu = P[3]
    A = P[4]
    sig = P[5]
    p = fs * p_peak(x, mu,  A, sig) + (1.-fs)*p_bkg(x, x0, x1, y0, y1)
    return p



# likelihood function for signal and background
def LH(parameters, x, x0, x1):
    # calculate likelihood, seep Barlow p.93 for Least Sq.
    # calculate norm
    fs = parameters[0]
    y0 = parameters[1]
    y1 = parameters[2]
    mu = parameters[3]
    A = parameters[4]
    sig = parameters[5]
    p = fs*p_peak(x, mu,   sig, A) + (1.-fs)*p_bkg(x, x0, x1, y0, y1)
    lgl = -np.sum(np.log(p))
    return lgl








#%%


def analyze_events( fname, Nf = 500, dy = 0.2, Nevents = -1, check = 0):
    # Nf number of events to be fitted per event
    # dy = fractional width in y
    # select range in y
    # Nevents number of events to analyze
    x0 =  0.90#0.8974868#0.8666207#0.8513027 #0.9
    x1 =  1.02#1.0209512 #1.0518173#1.0698533 #1.2
    y1_range = 2.
    #y2_range = 2*180.
    #x_range = 0.1997
    i = 1000
    npzfile = np.load(fname)
    xr = npzfile['metap'][:i]
    yr1 = npzfile['cost_etap'][:i]
    #yr2 = npzfile['phi'][:]
    
    #AccWeights = npzfile['AccWeights']
    #xr = pi0mass
    #yr1 = pi0costheta_GJ
    #yr2 = pi0phi_GJ
    
    xsave=[]
    FS=[]
    Y0=[]
    Y1=[]
    Mu = []
    #a = []
    #d = []
    #Sig = []
    #Sig1 = []
    # select peak region
    sel_p = (x0 < xr) & (xr < x1)
    xsel = xr[sel_p]
    print('xsel=',xsel)
    ysel1 = yr1[sel_p]
    #ysel2 = yr2[sel_p]
    
    qf = np.ones_like(xsel)
    fit_res= []
    if Nevents <= 0:
        Nevents = len(xsel) 
    for i,x in enumerate(xsel[:Nevents]):
        if (i%10000 == 0):
           print (i, " events fitted")
        y11 = ysel1[i]
        #y22 = ysel2[i]
        
        # find events in y neighbor hood
        sel_n = np.sqrt(((y11 - ysel1)/y1_range)**2)  <= dy
        #sel_n = np.sqrt(((y11 - ysel1)/y1_range)**2 +  ((y22 - ysel2)/y2_range)**2) <= dy
        # perform a fit to these events
        xf = xsel[sel_n][:Nf]
        
        
        
        #print("xf",xf)
        # perform a likelihood fit to the events and separate them
        
        fs = 0.5
        y0 = .0224
        y1 = -0.1139
        mu = 0.958
        sig = 0.0099
        A = 50.
        
        # constraints on the fit parameters
        #limits = [[0.,1.], [0.,1.e8], [0.,1.e8], [0.52, 0.56], [0.007, 0.008], [0.015, 0.025]]
        limits = [[0.,1.], [-1e8,1e8], [-1e8,1e8], [ 0.94, 0.98],  [0.0085,0.020],  [0.1,1e8]]
        #limits = [[0.,1.], [0.,1.e8], [0.,1.e8], [0.,1.e8], [0.,1.e8],  [0.015,0.025], [0.005,0.015]]
        # store parameters in an array
        fit_par = np.array([fs,y0,y1, mu, sig, A])
        # maximize likelihood
        fit_model = minimize(LH, fit_par, args = (xf, x0, x1),  bounds = limits, method='L-BFGS-B')
        #fit_model.
        #print (fit_model.chisqr)
        #fit_model = minimize(LH, fit_par, args = (xf),  bounds = limits, method='SLSQP')
        
        
        fs = fit_model['x'][0]
        y0 = fit_model['x'][1]
        y1 = fit_model['x'][2]
        mu = fit_model['x'][3]
        sig = fit_model['x'][4]
        A = fit_model['x'][5]
        
        #print('fs=',fs,'y0=',y0,"y1=",y1, "mu= ", mu, "A= ",A, "D= ",D, "sig=", sig, "sig1=",  sig1)
        
        #plt.plot(xf, p_signal(xf, (fs, y0,y1), x0, x1), label="fit");
        
        
        fit_res.append(fit_model)
        if not fit_model['success']:
            print ("bad fit for event : ", i)
        # now calculate weights
        ws = fs*p_peak(x,mu,  sig, A)
        wb = (1.-fs)*p_bkg(x, x0, x1, y0, y1)
        # calculate q-factor  
        q = ws/(ws+wb)
        qf[i] = q
        #print("q-value",q)
        
        
    
        #xsave.append(xf)
        FS.append(fs)
        Y0.append(y0)
        Y1.append(y1)
        Mu.append(mu)
        #Sig.append(sig)
        #Sig1.append(sig1)
    return [FS,Y0,Y1, Mu,   xsel, ysel1, qf, fit_res]

#%%
#FS,Y0,Y1, Mu,  xsel, ysel1, ysel2, qf, fit_res = analyze_events('/w/halld-scifs17exp/halld2/home/rupesh/halld/bkg_study/sep_qf_etap.npz', Nf = 1000)
#FS,Y0,Y1,Mu,  xsel, ysel1, ysel2, qf, fit_res = analyze_events('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/sep_qf_etap.npz', Nf = 10)
#FS,Y0,Y1,Mu,  xsel, ysel1,  qf, fit_res = analyze_events('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfac_analysis.npz', Nf = 500)

FS,Y0,Y1,Mu,  xsel, ysel1,  qf, fit_res = analyze_events('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18.npz', Nf = 100)

#FS,Y0,Y1,Mu,xf, xsave, xsel, ysel1, ysel2, qf, fit_res = analyze_events('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_mc/etaprimeqfactoramc.npz', Nf = 1000)


#%%

h = B.histo(xsel, bins = 25)
hs = B.histo(xsel, weights = qf, bins = 25)
hb = B.histo(xsel, weights = 1-qf, bins = 25)

















