#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 12:39:00 2021

@author: rupeshdotel
"""

import numpy as np
import LT.box as B
import matplotlib.pyplot as plt

#%%


#define function for momentum transfer t between beam and recoil system as proton
def t( thr_egy, max_egy, theta_X, n):
    """
    #e_beam = beam energey
    #e_x = resonance energy, m_X = mass of resonance
    #theta_x = angle between resonance and beam 
    #p_X = momentum of resonace 
    """
    e_X = np.linspace(thr_egy, max_egy, n)
    e_beam = np.linspace(8.2, 8.8, n)
    m_X = 1.6 # mass of resonace in GeV
   
    p_X = np.sqrt(e_X**2 - m_X**2)
    beta_X = p_X/e_X
    T = m_X**2 - 2* e_beam * e_X * (1 - beta_X * np.cos(np.radians(theta_X)))
    return -T, e_X, beta_X, p_X, e_beam 
    
#%%

mom_t, e_X, beta_X, p_X,  e_beam = t(1.6, 10, 5, 500)
sel_t = (0.1 < mom_t) & (mom_t < 0.7)
mom_t = mom_t[sel_t]
e_X = e_X[sel_t]
p_X = p_X[sel_t]

#%%
B.plot_exp(mom_t, e_X, plot_title = "resonance energy vs momemtum transfer t")
plt.figure();
B.plot_exp(mom_t, p_X, plot_title = "resonace momentum vs momemtum transfer t")

