#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 11:20:31 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B

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
h = B.histo(metap, bins = 40)
h.set_fit_list(fit = [ 'A', 'mean', 'sigma', 'b0', 'b1'])
fit = h.fit()
h.plot_exp()
h.plot_fit()

#%%

#try likelihood fit
#define a  gaussian peak 
def peak(x, A):
    #A = B.Parameter(700., 'A')
    x0 = B.Parameter(0.956, 'x0')
    sigma = B.Parameter(.014, 'sigma')
    y = A * np.exp(-(x - x0.value )**2/(2.*sigma.value**2))
    return y

#linear bkg
def bkg(x, b0, b1):
    y = b0 + b1 * x
    return y
    
#total function
def signal(x, A, b0, b1):
    return peak(x, A) + bkg(x, b0, b1)

#define the likelihood that will be minimized 
def lk(x, p):
    A = p[0]
    b0 = p[1]
    b1 = p[2]
    y = signal(x, A, b0, b1)
    lgl = -np.sum(np.log(y))
    return lgl      


#guess the initial parameters
A = 710.
b0 = -1240.
b1 = 1461.          
guess = np.array([A, b0, b1])
#limits = [[0.,1000.], [-2000, 2000], [0., 2000] ]

results = minimize(lk, guess, args = (metap)) #,  bounds = limits, options={'disp': True})

print(results)

A = results.x[0]
b0 = results.x[1]
b1 = results.x[2]

h.plot_exp()
B.plot_exp(metap, peak(metap, A))
B.plot_exp(metap, bkg(metap, b0, b1))
B.plot_exp(metap, signal(metap, A, b0, b1))



