#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 10:06:13 2021

@author: rupeshdotel

for testing only
"""

import sympy
from sympy import symbols, diff
import numpy as np

#%%
#define errors 
dA = 3.5
dm = 0.003
ds = 0.02
da0 = 0.004
da1 = 0.005 

#%%
#get the expression for 2 factor

x, A, m,  s, a0, a1  = symbols('x A m s a0 a1')

gaus_sig = A*sympy.exp(-(x - m)**2/(2.*s**2)) # signal
lin_bkg = a0 + a1*x #background

Total = gaus_sig + lin_bkg #total

q  = gaus_sig/(Total) #qfactor

#expression for partial derivatives with respect to different variables
dq_dA = diff(q, A)
dq_dm = diff(q, m)
dq_ds = diff(q, s)
dq_da0 = diff(q, a0)
dq_da1 = diff(q, a1)

#resultant expression for error on q-factors 
dq = sympy.sqrt((dq_dA * dA)**2 + (dq_dm * dm)**2 + (dq_ds * ds)**2 + (dq_da0 * da0)**2 +  (dq_da1 * da1)**2)




