#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 11:20:31 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import class_fit as cf

from scipy.optimize import minimize
import scipy.stats as stats

#%%
"""
N = 25
x = np.linspace(0, 20 , N)
e = np.random.normal(loc = 0.0, scale = 5.0, size = N) #introduce noise
slope = 3
y = np.abs( slope * x + e)

h = B.histo(bin_center = x, bin_content = y)

h.set_fit_list(['b0', 'b1'])
fit = B.linefit(x, y, e)

B.plot_exp(x, y, e)
B.plot_line(fit.xpl, fit.ypl)



#%%

def MLERegression(params):
    N = 25
    x = np.linspace(10, 20 , N)
    e = np.random.normal(loc = 0.0, scale = 5.0, size = N)
    y = np.abs(3 * x + e)
    intercept, slope, sd = params[0], params[1], params[2] # inputs are guesses at our parameters
    yhat = intercept + slope * x # predictions# next, we flip the Bayesian question
    # compute PDF of observed values normally distributed around mean (yhat)
    # with a standard deviation of sd
    negLL = -np.sum(stats.norm.logpdf(y, loc = yhat, scale = sd))# return negative LL
    return(negLL)

guess = np.array([0.1, 3, 5])

results = minimize(MLERegression, guess, method = 'Nelder-Mead', 
 options={'disp': True})

print(results)
offset = results.x[0]
slope = results.x[1]
err = results.x[2]

def func(X, M,  C):
    return M * X + C

B.plot_line(x, func(x, slope, offset))

"""

#%%
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events.npz')
metap = d['metap']


#chi square fit
h = B.histo(metap, bins = 20)
h.set_fit_list(fit = [ 'A', 'mean', 'sigma', 'b0', 'b1'])
fit = h.fit()
h.plot_exp()
h.plot_fit()


#%%

#try likelihood fit
#define a  gaussian peak 
"""
def peak(x):
    Aa = B.Parameter(1500., 'Aa')
    x0 = B.Parameter(0.957, 'x0')
    sigma = B.Parameter(.01665, 'sigma')
    #y = 1./(np.sqrt(2.*np.pi)*sigma.value) * np.exp(-(x - x0.value )**2/(2.*sigma.value**2))
    y = Aa.value * np.exp(-(x - x0.value )**2/(2.*sigma.value**2))
    return y

#linear bkg
def bkg(x, b0, c0):
    #y = b0 + b1 * x
    f = cf.gauss_fit()
         
    # set initial parameters for fitting
    #f.A.set(1500) 
    #f.x0.set(0.956)
    #f.sigma.set(0.01)
    f.b0.set(0.05)
    f.db0.set(0.5)
    f.c0.set(500)
    
    # set bounds for the fit parameters
    #f.A_min.set(0.); f.A_max.set(1e5)
    #f.x0_min.set(0.95); f.x0_max.set(0.96)
    #f.sigma_min.set(0.008); f.sigma_max.set(0.016)
    f.c0_min.set(0.); f.c0_max.set(1e5)
    f.b0_min.set(0.00); f.b0_max.set(0.20)
    y = f.bt_bkg(x)
    return y
    
#total function
def signal(x, b0, c0):
    return peak(x) + bkg(x, b0, c0)

"""

f = cf.gauss_fit()
         
# set initial parameters for fitting
f.A.set(1500.) 
f.x0.set(0.956)
f.sigma.set(0.01)
f.b0.set(0.05)
f.db0.set(0.5)
f.c0.set(500)

# set bounds for the fit parameters
f.A_min.set(0.); f.A_max.set(1e5)
f.x0_min.set(0.95); f.x0_max.set(0.96)
f.sigma_min.set(0.008); f.sigma_max.set(0.016)
#f.c0_min.set(0.); f.c0_max.set(1e5)
f.b0_min.set(0.00); f.b0_max.set(0.20)


#define the likelihood that will be minimized 
def lk(x, p):
    #A = p[0]
    #K = np.zeros_like(x)
    b0 = p[0]
    c0 = p[1]
    #b0 = f.b0.value
    #c0 = f.c0.value
    
    
    y = f.signal_bt_bkg(x)
    #for i in p:
        #if i > 0:
           # k = np.log(y)
    lgl = -np.sum(np.log(y))
    return lgl      


#guess the initial parameters
#A = 1500.
#b0 = -8000.
#b1 = 9500.          
#guess = np.array([ b0, b1])
#limits = [[-1e5, 1e5], [-5e5, 5e5]]
#q = cf.gauss_fit()
#b0 = q.b0.set(0.05)
#c0 = q.c0.set(500)

b0 = f.b0.value
c0 = f.c0.value


fit_par = np.array([b0, c0])
results = minimize(lk, fit_par, args = (metap),   method='L-BFGS-B')

print(results)

#%%
#A = results.x[0]
b0 = results.x[0]
b1 = results.x[1]


h.plot_exp()

B.plot_exp(metap, peak(metap))
B.plot_exp(metap, bkg(metap, b0, b1))
B.plot_exp(metap, signal(metap, b0, b1))



