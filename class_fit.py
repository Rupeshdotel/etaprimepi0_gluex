#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:34:13 2020

@author: rupeshdotel
"""

from mpl_toolkits import mplot3d
import numpy as np
import LT.box as B
import matplotlib.pyplot as plt



#%%

class gauss_bt_fit():
    
     def __init__(self, A = 10000., x0 = 0.956, sigma = 0.01,  c0 = 5e4, 
                  
                  b0 = 0.14, db0 = 0.50, m0 = 0.86, m1 = 1.025):
        # peak parameters
        self.A = B.Parameter(A, 'A')
        self.x0 = B.Parameter(x0, 'x0')
        self.sigma = B.Parameter(sigma, 'sigma')
        # back ground parameters
        
        self.c0 = B.Parameter(c0, 'c0')
        
        self.b0 = B.Parameter(b0, 'b0')
        self.db0 = B.Parameter(db0, 'db0')
        self.m0 = B.Parameter(m0, 'm0')
        self.m1 = B.Parameter(m1, 'm1')
        #self.b1 = B.Parameter(b1, 'b1')
        # fit list, which parameters are to be varied
        self.fit_par = {
        "A" : self.A, 
        "sigma":  self.sigma, \
        "c0" : self.c0, \
        "b0" : self.b0 }
        
        self.set_fit_list()
        self.fit_list = [self.A,  self.sigma,  self.c0, self.b0]
                         #self.b0, self.db0,
                         #self.b2, self.db2]
        
     def gauss(self, x):
       
        y = self.A()*np.exp(-((x-self.x0())/(2.*self.sigma()))**2)
        return y
     
    
     def B23(self, x):
        y = 3*x*x*(1-x)
        return y
     
     def B34(self, x):
        y = 4*x*x*x*(1-x)
        return y
    
     
     
     
     def b34(self, c0, x, m0, m1):
        self.b1 = self.b0() + self.db0()
        self.a = self.db0()/(self.m1() - self.m0())
        self.b = (self.b1 * self.m0() - self.b0() * self.m1())/(self.m0() - self.m1())
        self.X = self.a * x + self.b
        y = 4 * self.X**3 * (1 - self.X)
        return self.c0() * y
    
    
     
     def bkg(self, x):
        #   bernstein polynomial background
        y =   self.b34(self.c0, x, self.m0, self.m1) 
        return y
    
     
     def signal(self,x):
        
        return self.gauss(x) + self.bkg(x)

     def fit(self, x, y, dy):
        self.x = x
        self.y = y
        if dy is None:
            self.dy = np.ones_like(y)
        else:
            self.dy = dy
        if len(self.fit_list) == 3:
            self.pl = np.array([0.   ,     0.009,       0. ])
            self.pu = np.array([1e5 ,      0.016,     1e5])
        if len(self.fit_list) == 4:
            self.pl = np.array([ 0.   ,     0.008,      0., 0.0 ])
            self.pu = np.array([ 1e5 ,      0.016,    1e5, 0.2 ])
        self.fit_res = B.genfit(self.signal, self.fit_list, self.x, self.y,
                                y_err = self.dy, bounds = (self.pl, self.pu))
        
     def plot_fit(self):
        B.plot_line(self.fit_res.xpl, self.fit_res.ypl)
        
     #function to get the correlation matrix elements i and j are the parameter indices
     def calc_corr(self, i, j):
        self.i = i
        self.j = j
        #if i == j:
         #   return self.fit_res.covar[i][i]/self.fit_res.parameters[i].err**2
        #else:
        return self.fit_res.covar[i][j]/(self.fit_res.parameters[i].err*self.fit_res.parameters[j].err)
    
     
     def get_corr(self):
        self.z = self.fit_res.covar
        self.rho = np.zeros_like(self.z)
        
        for i in range(self.z.shape[0]):
            for j in range(self.z.shape[1]):
                self.rho[i,j] = self.z[i, j]/np.sqrt(self.z[i, i] * self.z[j, j])
        
         
        
        return self.rho
        
    
     def show_fit_list(self):
        """
        Show the current fit list

        Returns
        -------
        None.

        """
        print("\nCurrent fit list : ", [k.name for k in self.fit_list])
        print("\nAvailable parameters: [ 'A',  'sigma',  'c0', 'b0']")
        
        
     def set_fit_list(self, fit = [ 'A',  'sigma'] ):
        """

        Define which parameters are to be fitted.

        The default list is ::
        
           fit = [ 'A', 'sigma']

        to use all parameters::

           h.set_fit_list( fit = [ 'A', 'sigma', 'c0', 'b0'])
           
        """
        if fit==[]:
            # empty list show possibilities
            print('possible fit parameters:')
            print(list(self.fit_par.keys()))
            return
        self.fit_names = fit
        self.fit_list = []
        for key in self.fit_names:
            try:
                curr_par_name = self.fit_par[key]
            except:
                print('cannot use parameter :', key, ' (does not exist ?!)')
                continue
            self.fit_list.append(curr_par_name)
        # end of fitting list
        