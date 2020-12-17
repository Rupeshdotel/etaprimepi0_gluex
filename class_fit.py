#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:34:13 2020

@author: rupeshdotel
"""


import numpy as np
import LT.box as B


#%%

class gauss_fit():
    
     def __init__(self, A_min = 10., A = 10000.,  A_max = 1e9, 
                        x0_min = 0.94,    x0 = 0.956, x0_max = 0.98,
                        sigma_min = 0.005,   sigma = 0.01, sigma_max = 0.30,
                        c0_min = 10,  c0 = 1e4, c0_max = 5e4, 
                        b0_min = 10,  b0 = 0.14, b0_max = 0.20, 
                        db0 = 0.50, 
                        k0_min = 10,  k0 = 100., k0_max = 500., 
                        k1_min = 10,  k1 = 100., k1_max = 500., 
                        m0 = 0.86, m1 = 1.04):
        # peak parameters
        self.A_min = B.Parameter(A_min, 'A_min')
        self.A = B.Parameter(A, 'A')
        self.A_max = B.Parameter(A_max, 'A_max')
        
        self.x0_min = B.Parameter(x0_min, 'x0_min')
        self.x0 = B.Parameter(x0, 'x0')
        self.x0_max = B.Parameter(x0_max, 'x0_max')
        
        self.sigma_min = B.Parameter(sigma_min, 'sigma_min')
        self.sigma = B.Parameter(sigma, 'sigma')
        self.sigma_max = B.Parameter(sigma_max, 'sigma_max')
        
        # back ground parameters
        
        self.c0_min = B.Parameter(c0_min, 'c0_min')
        self.c0 = B.Parameter(c0, 'c0')
        self.c0_max = B.Parameter(c0_max, 'c0_max')
        
        self.b0_min = B.Parameter(b0_min, 'b0_min')
        self.b0 = B.Parameter(b0, 'b0')
        self.b0_max = B.Parameter(b0_max, 'b0_max')
        
        self.db0 = B.Parameter(db0, 'db0')
        
        self.k0_min = B.Parameter(k0_min, 'k0_min')
        self.k0 = B.Parameter(k0, 'k0')
        self.k0_max = B.Parameter(k0_max, 'k0_max')
        
        self.k1_min = B.Parameter(k1_min, 'k1_min')
        self.k1 = B.Parameter(k1, 'k1')
        self.k1_max = B.Parameter(k1_max, 'k1_max')
        
        self.m0 = B.Parameter(m0, 'm0')
        self.m1 = B.Parameter(m1, 'm1')
        
        
        #self.b1 = B.Parameter(b1, 'b1')
        # fit list, which parameters are to be varied
        self.fit_par = {
        "A" : self.A, 
        "x0" : self.x0, \
        "sigma":  self.sigma, \
        "c0" : self.c0, \
        "b0" : self.b0, \
        "k0" : self.k0, \
        "k1" : self.k1 }
        
        self.set_fit_list()
        self.fit_list = [self.A, self.x0,  self.sigma,  self.c0, self.b0, self.k0, self.k1 ]
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
    
     
     
     
     def bt_bkg(self,  x):
        self.b1 = self.b0() + self.db0()
        self.a = self.db0()/(self.m1() - self.m0())
        self.b = (self.b1 * self.m0() - self.b0() * self.m1())/(self.m0() - self.m1())
        self.X = self.a * x + self.b
        y = 4 * self.X**3 * (1 - self.X)
        return self.c0() * y
    
    
     
     def lin_bkg(self, x):
        
        y = self.k0() + self.k1()*x 

        return y
    
     
     def signal_bt_bkg(self, x):
        
        return self.gauss(x) + self.bt_bkg(x)
    
     def signal_lin_bkg(self,x):
        
        return self.gauss(x) + self.lin_bkg(x)

     def fit_gaussbt(self, x, y, dy):
        self.x = x
        self.y = y
        if dy is None:
            self.dy = np.ones_like(y)
        else:
            self.dy = dy
        
        
        self.pl = np.array([ self.A_min() , self.x0_min(),  self.sigma_min(), self.c0_min(), self.b0_min() ])
        self.pu = np.array([ self.A_max() , self.x0_max(),  self.sigma_max(), self.c0_max(), self.b0_max() ])
        
        self.fit_res = B.genfit(self.signal_bt_bkg, self.fit_list, self.x, self.y,
                                y_err = self.dy, bounds = (self.pl, self.pu))
        
     def fit_gausslin(self, x, y, dy):
        self.x = x
        self.y = y
        if dy is None:
            self.dy = np.ones_like(y)
        else:
            self.dy = dy
        
        
        
        self.pl = np.array([ self.A_min() , self.x0_min(),  self.sigma_min(), self.k0_min(), self.k1_min()  ])
        self.pu = np.array([ self.A_max() , self.x0_max(),  self.sigma_max(), self.k0_max(), self.k1_max() ])
        self.fit_res = B.genfit(self.signal_lin_bkg, self.fit_list, self.x, self.y,
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
        print("\nAvailable parameters: [ 'A', 'x0',  'sigma',  'c0', 'b0', 'k0', 'k1']")
        
        
     def set_fit_list(self, fit = [ 'A',  'sigma'] ):
        """

        Define which parameters are to be fitted.

        The default list is ::
        
           fit = [ 'A', 'sigma']

        to use bernstien bkg  use parameters::

           h.set_fit_list( fit = [ 'A', 'x0',  'sigma', 'c0', 'b0'])
           
        to use linear bkg  use parameters::

           h.set_fit_list( fit = [ 'A', 'x0', 'sigma', 'k0', 'k1'])
           
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
#end of class

#%%
#if __name__ == "__main__":
 #  gauss_bt_fit()



        