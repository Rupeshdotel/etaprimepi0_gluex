#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 20:11:52 2021

@author: rupeshdotel
"""

#qfactors with least square fitting

import numpy as np
import matplotlib.pyplot as plt
import LT.box as B
import scipy.spatial as SP
from root_numpy import array2tree, array2root, tree2array
import ROOT as R
from root_numpy import tree2array

#%%

A = B.Parameter(50., 'A')
x0 = B.Parameter(0.956, 'x0')
sigma = B.Parameter(.005, 'sigma')

def gaus(x):
    return A()*np.exp(-(x - x0() )**2/(2.*sigma()**2))


a0 = B.Parameter(1., 'a0')
a1 = B.Parameter(1., 'a1')


def lin_bkg(x):
    return a0() + a1()*x 

def signal(x):
    return gaus(x) + lin_bkg(x)

#%%

#f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_unique_cons_17_18_incamo.npz')

f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_unique_cons.npz')
metappi0 = f['metappi0']

#%%
# cos_gj is the costheta in Gottfried Jackson frame
cos_theta_etap_gj = f['cost_etap_gj']


#%%
phi_etap_gj = f['phi_etap_gj'] #phi is in degrees
#get cosphi and sinphi 

#np.cos(theta), numpy requires theta in radians, np.radians(theta) converts degrees to radians
cos_phi_etap_gj = np.cos(np.radians(phi_etap_gj)) 
sin_phi_etap_gj = np.sin(np.radians(phi_etap_gj)) 


#%%
# m_etap is the mass of etaprime  
m_etap = f['metap']


#%%
# M_min and M_max are the left and right edge values
M_min = 0.90; M_max = 1.02
sel = (M_min < m_etap) & (m_etap < M_max)
# M_inv_sel are the selected events within the range
M_inv_sel = m_etap[sel]

#%%
#k = 100
#pair costheta with itself for distance calculation
cos_theta_etap_gj_sel_i = cos_theta_etap_gj[sel][:]
cos_theta_etap_gj_sel_j = cos_theta_etap_gj[sel][:]

#pair cosphi with sinphi for distance calculation
cos_phi_etap_gj_sel_i = cos_phi_etap_gj[sel]
sin_phi_etap_gj_sel_j = sin_phi_etap_gj[sel]



cosphi_sinphi_etap_gj_pair = np.array([cos_phi_etap_gj_sel_i ,sin_phi_etap_gj_sel_j]).T 


# array of angle pairs
cos_theta_etap_gj_pair = np.array([cos_theta_etap_gj_sel_i ,cos_theta_etap_gj_sel_j]).T 


# calculate the distance between all gj angles in array 
dcos_theta_etap_gj_a = (SP.distance.pdist(cos_theta_etap_gj_pair))/2 #divide by 2 to correct for double counting

#calculate the distance between all gj phi angles in array 
dphi_etap_gj_a = SP.distance.pdist(cosphi_sinphi_etap_gj_pair)

# cos_theta_gj_range is the total range of costheta variable
cos_theta_gj_range = 4.0

# phi_gj_range is the maximum distance possible between (cosphi, sinphi) points 
phi_gj_range = 8.0

# convert distance array  to symmetric matrix and normalize
dcos_theta_etap_gj_m = SP.distance.squareform(dcos_theta_etap_gj_a)/cos_theta_gj_range
dphi_etap_gj_m = SP.distance.squareform(dphi_etap_gj_a)/phi_gj_range






#normalize with maximum distance between variables
#dcos_theta_etap_gj_m_norm = dcos_theta_etap_gj_m
#dphi_etap_gj_m_norm = dphi_etap_gj_m

#%%
# Nf are  the number of neighboring events we choose
n_near = 200


# k is an index for testing
#k = 250

qf = np.ones_like(M_inv_sel[:])

# X0 are the mean values from the fit
X0 = np.ones_like(M_inv_sel[:])

# S0 are the sigma values from the fit
S0 = np.ones_like(M_inv_sel[:])
do_plot  = False


q_err = np.ones_like(M_inv_sel[:])


#%%
for i,M_inv_loc in enumerate(M_inv_sel[2561:2563]):
    
    # select neghboring events for the current event combine distance of costheta and phi normalized by  pythagorian theorem   
    i_near = np.argsort(np.sqrt((dcos_theta_etap_gj_m[i]**2 +  dphi_etap_gj_m[i]**2 )))[:n_near]
    
   
        
    M_inv_neighbor = M_inv_sel[i_near]
   
    
    h = B.histo(M_inv_neighbor, bins = 15)
    h_sel = h.bin_content > 0
    M = h.bin_center[h_sel]
    C = h.bin_content[h_sel]
    dC = h.bin_error[h_sel]
    
    A = B.Parameter(C.max(), 'A')
    x0 = B.Parameter(0.956, 'x0')
    sigma = B.Parameter(.005, 'sigma')
    a0 = B.Parameter(1., 'a0')
    a1 = B.Parameter(1., 'a1')
    
    fit = B.genfit(signal, [A, x0, sigma,  a0, a1],
                   x = M, y = C, y_err = dC, print_results = False )
    
    if do_plot:
        plt.figure()
        h.plot_exp()
        B.plot_line(fit.xpl, fit.ypl)
        mr = np.linspace(M[0], M[-1], 1000)
        B.plot_line(mr, gaus(mr))
        B.plot_line(mr, lin_bkg(mr))
    
    A = B.Parameter(fit.parameters[0].value, 'A') #amplitude value
    dA = B.Parameter(fit.parameters[0].err, 'dA') #ampltiude error
    
    x0 = B.Parameter(fit.parameters[1].value, 'x0') #mean value
    dx0 = B.Parameter(fit.parameters[1].err, 'dx0') #mean error
    
    sigma = B.Parameter(fit.parameters[2].value, 'sigma') # sigma value
    dsigma = B.Parameter(fit.parameters[2].err, 'dsigma') # sigma error
    
    a0 = B.Parameter(fit.parameters[3].value, 'a0') # a0 value
    da0 = B.Parameter(fit.parameters[3].err, 'da0') # a0 error
    
    a1 = B.Parameter(fit.parameters[4].value, 'a1') # a1 value
    da1 = B.Parameter(fit.parameters[4].err, 'da1') # a1 error
    
    
    ws = gaus(M_inv_loc)
    wb = lin_bkg(M_inv_loc)
    # calculate q-factor  
    q_loc = ws/(ws+wb)
    qf[i] = q_loc
    #X0[i] = x0.value
    #S0[i] = sigma.value
    
    #calculate errors on q-factors
    p = ((M_inv_loc - x0.value)/(sigma.value))**2 #power expression
    t = (A.value * np.exp(-0.5 * p)) + a0.value + a1.value # total function in denmoinator
    
    
    #these expressions can be calculated by hand or use sympy package in a different script
    dq_dA = (-A.value * np.exp(-p))/t**2 + (np.exp(-p/2)/t)
    
    dq_dx0 = ((A.value**2 * (x0.value - M_inv_loc)*np.exp(-p))/(sigma.value**2 * t**2)) - \
    (A.value * (x0.value - M_inv_loc)*np.exp(-p/2)/(sigma.value**2 * t))
    
    dq_dsigma = ((-A.value**2 * (M_inv_loc - x0.value)**2 *np.exp(-p))/(sigma.value**3 * t**2)) + \
    ((A.value * (M_inv_loc - x0.value)**2 *np.exp(-p/2))/(sigma.value**3 * t))
    
    dq_da0 = (-A.value * np.exp(-p/2))/(t**2)
    dq_da1 = (-A.value * M_inv_loc * np.exp(-p/2))/(t**2)
    
    #get covariance matrix elements
    cov_A_x0 = fit.covar[0][1]
    cov_A_sigma = fit.covar[0][2]
    cov_A_a0 = fit.covar[0][3]
    cov_A_a1 = fit.covar[0][4]
    
    cov_x0_sigma = fit.covar[1][2]
    cov_x0_a0 = fit.covar[1][3]
    cov_x0_a1 = fit.covar[1][4]
    
    cov_sigma_a0 = fit.covar[2][3]
    cov_sigma_a1 = fit.covar[2][4]
    
    cov_a0_a1 = fit.covar[3][4]
    
    
    
    
    
    #get the error
    dq_loc = np.sqrt(
        (dq_dA * dA.value)**2 + (dq_dx0 * dx0.value)**2 + (dq_dsigma * dsigma.value)**2 + \
    (dq_da0 * da0.value)**2 + (dq_da1 * da1.value)**2 \
        
     + 2* ( dq_dA * dq_dx0 * cov_A_x0 #cross terms
     +  dq_dA * dq_dsigma * cov_A_sigma
     +  dq_dA * dq_da0 * cov_A_a0 
     +  dq_dA * dq_da1 * cov_A_a1 
         
       +  dq_dx0 * dq_dsigma * cov_x0_sigma
        +   dq_dx0 * dq_da0 * cov_x0_a0
         +   dq_dx0 * dq_da1 * cov_x0_a1
          
          
           +   dq_dsigma * dq_da0 * cov_sigma_a0 
           +   dq_dsigma * dq_da1 * cov_sigma_a1
           
           +   dq_da0 * dq_da1 * cov_a0_a1 ) ) #uncertainty in a function of several variables 
                                                             #page 73 John R. Taylor Introduction to error analysis
                
    q_err[i] = dq_loc
    
    
    #get correlation matrix
    z = fit.covar
    corr_mtx = np.zeros_like(z) 
     
    
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
                corr_mtx[i,j] = z[i, j]/np.sqrt(z[i, i] * z[j, j])
   
        
    #check bad fits if any
    if((q_loc < 0) | (q_loc >1)):
        print(f'badfits : qfactor =  {q_loc:.3f}, q_error = {dq_loc:.3f}, m = {M_inv_loc:.3f},\
              ws = {gaus(M_inv_loc):.3f}, wb = {lin_bkg(M_inv_loc):.3f}, \
               cm = {corr_mtx}')
        plt.figure()
        h.plot_exp()
        B.plot_line(fit.xpl, fit.ypl)
        mr = np.linspace(M[0], M[-1], 1000)
        B.plot_line(mr, gaus(mr))
        B.plot_line(mr, lin_bkg(mr))
        fit = B.genfit(signal, [A, x0, sigma,  a0, a1],
                   x = M, y = C, y_err = dC, print_results = True )
        plt.title("") 
        plt.xlabel("")
        plt.ylabel("")
        
        
        
#%%

sel_q = (qf<1) & (qf>0)
qsel = qf[sel_q]
metap_min = 0.90
metap_max = 1.02
metap_binwidth = 3e-3 # bin width in GeV
bins_metap = int((metap_max - metap_min)/metap_binwidth)
h  = B.histo(M_inv_sel[sel_q], bins = bins_metap)   
hs  = B.histo(M_inv_sel[sel_q], weights = qsel,  bins = bins_metap)   
hb  = B.histo(M_inv_sel[sel_q], weights = 1-qsel,  bins = bins_metap)   


plt.figure()
h.plot_exp()
hs.plot_exp()
hb.plot_exp()
plt.title("Q FACTOR weighted histogram") 
plt.xlabel("")
plt.ylabel("")
    
#%%

cost_min = -1.0
cost_max = 1.0
bins_cost = 20

metappi0_min = 1.0
metappi0_max = 2.5
metappi0_binwidth = 50e-3 # bin width in GeV
bins_metappi0 = int((metappi0_max - metappi0_min)/metappi0_binwidth)

h2_ep_cost_weighted = B.histo2d(metappi0[sel][sel_q], cos_theta_etap_gj[sel][sel_q],  bins = (bins_metappi0, bins_cost), 
                            range = [[metappi0_min, metappi0_max],[cost_min, cost_max]],
                       title = "Angular distribution in GJ Frame",
               weights = qsel)


#%%

np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qweights.npz',
         qf = qf)

#%%
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qweights_2d.npz',
         qf = qf)

#test github desktop
#%%


px_pr = f['px_pr']
px_etapr = f['px_etapr']
px_pi0 = f['px_pi0']


py_pr = f['py_pr']
py_etapr = f['py_etapr']
py_pi0 = f['py_pi0']


pz_pr = f['pz_pr']
pz_etapr = f['pz_etapr']
pz_pi0 = f['pz_pi0']

e_pr = f['e_pr']
e_etapr = f['e_etapr']
e_pi0 = f['e_pi0']



px_beam = f['px_beam']
py_beam = f['py_beam']
pz_beam = f['pz_beam']
e_beam = f['e_beam']


pol = f['pol']
#%%



px_pr = px_pr[sel][sel_q]
px_etapr = px_etapr[sel][sel_q]
px_pi0 = px_pi0[sel][sel_q]

py_pr = py_pr[sel][sel_q]
py_etapr = py_etapr[sel][sel_q]
py_pi0 = py_pi0[sel][sel_q]


pz_pr = pz_pr[sel][sel_q]
pz_etapr = pz_etapr[sel][sel_q]
pz_pi0 = pz_pi0[sel][sel_q]

e_pr = e_pr[sel][sel_q]
e_etapr = e_etapr[sel][sel_q]
e_pi0 = e_pi0[sel][sel_q]




px_beam = px_beam[sel][sel_q]
py_beam = py_beam[sel][sel_q]
pz_beam = pz_beam[sel][sel_q]
e_beam = e_beam[sel][sel_q]


pol = pol[sel][sel_q]


#%%

np.savez("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_qfactor.npz",
         
         
        px_pr = px_pr,
        px_etapr = px_etapr,
        px_pi0 = px_pi0,


        py_pr = py_pr,
        py_etapr = py_etapr,
        py_pi0 = py_pi0,


        pz_pr = pz_pr,
        pz_etapr = pz_etapr,
        pz_pi0 = pz_pi0,

        e_pr = e_pr,
        e_etapr = e_etapr,
        e_pi0 = e_pi0,



        px_beam = px_beam,
        py_beam = py_beam,
        pz_beam = pz_beam,
        e_beam = e_beam,


        pol = pol,
        qf = qsel
         )


#%%























"""
#for undergrad project only

px_pi0 = f['px_pi0']
py_pi0 = f['py_pi0']
pz_pi0 = f['pz_pi0']
e_pi0 = f['e_pi0']


px_etapr = f['px_etapr']
py_etapr = f['py_etapr']
pz_etapr = f['pz_etapr']
e_etapr = f['e_etapr']


px_pr = f['px_pr']
py_pr = f['py_pr']
pz_pr = f['pz_pr']
e_pr = f['e_pr']

#%%
px_pi0 = px_pi0[sel][sel_q]
py_pi0 = py_pi0[sel][sel_q]
pz_pi0 = pz_pi0[sel][sel_q]
e_pi0 = e_pi0[sel][sel_q]

px_etapr = px_etapr[sel][sel_q]
py_etapr = py_etapr[sel][sel_q]
pz_etapr = pz_etapr[sel][sel_q]
e_etapr = e_etapr[sel][sel_q]

px_pr = px_pr[sel][sel_q]
py_pr = py_pr[sel][sel_q]
pz_pr = pz_pr[sel][sel_q]
e_pr = e_pr[sel][sel_q]

cos_theta_etap_gj = cos_theta_etap_gj[sel][sel_q]
phi_etap_gj = phi_etap_gj[sel][sel_q]

#%%
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_four_vec.npz',
         
        px_pi0 = px_pi0,
        py_pi0 = py_pi0,
        pz_pi0 = pz_pi0,
        e_pi0 = e_pi0,

        px_etapr = px_etapr,
        py_etapr = py_etapr,
        pz_etapr = pz_etapr,
        e_etapr = e_etapr,

        px_pr = px_pr,
        py_pr = py_pr,
        pz_pr = pz_pr,
        e_pr = e_pr)



#%%


def npz_to_4vecs(filename):

    
    #intree = filename.Get('qfactortree')
    #f = tree2array(intree)
    f = filename
    
    px_pi0 = f['px_pi0']
    py_pi0 = f['py_pi0']
    pz_pi0 = f['pz_pi0']
    e_pi0 = f['e_pi0']
    
    
    px_etapr = f['px_etapr']
    py_etapr = f['py_etapr']
    pz_etapr = f['pz_etapr']
    e_etapr = f['e_etapr']
    
    
    px_pr = f['px_pr']
    py_pr = f['py_pr']
    pz_pr = f['pz_pr']
    e_pr = f['e_pr']
    
    '''
    px_beam = f ['px_beam']
    py_beam = f['py_beam']
    pz_beam = f['pz_beam']
    e_beam = f['e_beam']
    '''
    
    
    p4_pi0 = []
    p4_etaprime = []
    p4_proton = []
    
    
    
    for i in range(len(px_pr)):
        
        #pi0
        wx = px_pi0[i]
        wy = py_pi0[i]
        wz = pz_pi0[i]
        we = e_pi0[i]
        
        #etaprime
        vx = px_etapr[i]
        vy = py_etapr[i]
        vz = pz_etapr[i]
        ve = e_etapr[i]
        
        #proton
        ux = px_pr[i]
        uy = py_pr[i]
        uz = pz_pr[i]
        ue = e_pr[i]
        
        
        '''
        #4 vector combonents for beam
        be = e_beam[i]
        bx = px_beam[i]
        by = py_beam[i]
        bz = pz_beam[i]'''
    
        #x, y, z  and e of final states combined 
        pi0 = np.array([wx, wy, wz, we], dtype = 'float32')
        etapr = np.array([vx, vy, vz, ve], dtype = 'float32')
        proton = np.array([ux, uy, uz, ue], dtype = 'float32')
        
        
        
        
        
        p4_pi0.append(pi0)
        p4_etaprime.append(etapr)
        p4_proton.append(proton)
        
    return p4_pi0, p4_etaprime, p4_proton

#%%

fn  = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_four_vec.npz')

p4_pi0, p4_etaprime, p4_proton = npz_to_4vecs(fn)

#%%
p4_pi0 = np.array(p4_pi0)
p4_etaprime = np.array(p4_etaprime)
p4_proton = np.array(p4_proton)

#%%
p4_pi0.dtype = [('P4_Pi0', np.float32, (4,))]
p4_etaprime.dtype = [('P4_Etaprime', np.float32, (4,))]
p4_proton.dtype = [('P4_Proton', np.float32, (4,))]




#%%
qs_weights = qsel
qs_weights = qs_weights.astype([('q_factors', np.float32, (1,))])
cos_theta_etap_gj = cos_theta_etap_gj.astype([('Etaprime_costheta_gj', np.float32, (1,))])
phi_etap_gj = phi_etap_gj.astype([('Etaprime_phi_gj', np.float32, (1,))])


#%%



array2root(p4_pi0, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='finalstate_4vecs', mode = 'recreate')

array2root(p4_etaprime, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='finalstate_4vecs', mode = 'update')
array2root(p4_proton, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='finalstate_4vecs', mode = 'update')

array2root(cos_theta_etap_gj, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/gluexI_etaprimepi0_qfactors.root", treename='finalstate_4vecs', mode = 'update')

array2root(phi_etap_gj, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/gluexI_etaprimepi0_qfactors.root", treename='finalstate_4vecs', mode = 'update')

#%%


#test the file 
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/gluexI_etaprimepi0_qfactors.root")
intree = rfile.Get('finalstate_4vecs')
d = tree2array(intree)

#%%
etaprime = d['P4_Etaprime']

#get the component in arrays
etaprime_px = etaprime.T[0]
etaprime_py = etaprime.T[1]
etaprime_pz = etaprime.T[2]
etaprime_e = etaprime.T[3]
#%%

#get the invariant mass (etaprime)
m = np.sqrt(etaprime_e**2 - etaprime_px**2 - etaprime_py**2 - etaprime_pz**2)

etaprime_qfactors = d['q_factors']
etaprime_qfactors = etaprime_qfactors.flatten()

"""


 

