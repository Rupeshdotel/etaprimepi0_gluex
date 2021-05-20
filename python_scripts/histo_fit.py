#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 16:05:33 2020

@author: boeglinw
"""

import numpy as np
import LT.box as B
from scipy import integrate

#%%

# load the data
#d = B.get_file('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/pi0_accid.data')
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/fitpi0.npz')
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/pi0qfactor_0303.npz')

#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/etapdata.npz')
#mass = d['metap']

#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_18_08_pwa.npz')
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/cons_17_18_pol_0.npz')
mass = d['metap']

#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/pi0data.npz')
#mass = d['mpi0']


h = B.histo(mass, bins = 20, range = [0.85, 1.05])
#%%
#h = B.histo(mass, bins = 20)
M = h.bin_center
C = h.bin_content
dC = h.bin_error
dM = M[1] - M[0]

# get the bin width
#dM = d.par.get_value('dx')

# get the histogram data
#M = d['xb']
#C = d['cont']
#dC = d['dcont']
#dC = np.sqrt(C)

# plot the data

B.plot_exp(M, C, dC)




# setup
# simple gaussian
A = B.Parameter(C.max(), 'A')
#x0 = B.Parameter(0.135, 'x0')
#x0 = B.Parameter(0.135, 'x0')
x0 = B.Parameter(0.956, 'x0')
#x01 = B.Parameter(0.545, 'x01')
sig = B.Parameter(.005, 'sigma')


def gaus(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2))

# second gauss for a doublegauss
D = B.Parameter(100., 'D')
sig1 = B.Parameter(.005, 'sigma1')


#def doublegaus(x):
#    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2)) + D()*np.exp(-(x - x0() )**2/(2.*sig1()**2))



a0 = B.Parameter(1., 'a0')
a1 = B.Parameter(1., 'a1')
a2 = B.Parameter(1., 'a2')
#a2 = B.Parameter(207314., 'a2')



# background simple polynomial
def lin_bkg(x):
    #return a0() + a1()*x + a2()*x*x
    return a0() + a1()*x 


# fit function

def signal(x):
    return gaus(x) + lin_bkg(x)



'''
def gaus(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2))

Tau = B.Parameter(0.0225,'Tau')
def Voigt(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sig()**2)) + (D()* 0.25*Tau()**2/((x - x0())**2 + 0.25*Tau()**2))

def BW(x):
    return (D()* 0.25*Tau()**2/((x - x0())**2 + 0.25*Tau()**2))

a2 = B.Parameter(0., 'a2')
def pol_bkg(x):
    return a0() + a1()*x + a2()*x**2

alpha = B.Parameter(10000., 'alpha')
beta = B.Parameter(10000., 'beta')
gamma = B.Parameter(10., 'gamma')
delta = B.Parameter(1., 'delta')

def tanh_bkg(x):
    return alpha() + beta()*np.tanh(gamma()*(x - delta()))
 

'''

# select fit range

##%%

#M_min = 0.10; M_max = 0.18
M_min = 0.90; M_max = 1.02
sel = (M_min < M) & (M < M_max)


fit = B.genfit(signal, [A, x0, sig,  a0, a1],
               x = M[sel], y = C[sel], y_err = dC[sel] )

# plot the fit
B.plot_line(fit.xpl, fit.ypl)
B.pl.xlabel(r'$M_{\gamma\gamma}$')
B.pl.ylabel('Counts')
B.pl.title(' Gaussian signal with linear background')
mr = np.linspace(M[0], M[-1], 1000)
B.plot_line(mr, gaus(mr))

#%% number of counts by integrating the fit function

#l = 0.935
#u = 0.980

mean = fit.parameters[1].value
sigma = fit.parameters[2].value

#for intergration take mean +- 3 sigma
l = mean - 3 * sigma
u = mean + 3 * sigma

T_fit = integrate.quad(fit.func, l, u)[0]/dM
B_fit = integrate.quad(lin_bkg, l, u)[0]/dM
S_fit = integrate.quad(gaus, l, u)[0]/dM

print(f"\nTotal fitted counts = {T_fit:.0f},  Signal Counts = {S_fit:.0f}, Back Ground = {B_fit:.0f} ")

# total exp. counts

sel_data = (l < M) & (M < u)
C_tot = C[sel_data].sum()

print(f"\nExp. to Fit ratio = {C_tot/T_fit:.3f}")
