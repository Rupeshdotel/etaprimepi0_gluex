#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:32:17 2021

@author: rupeshdotel

generate kin tree for PWA from clean samples in the form of numpy arrays

"""

from root_numpy import array2tree, array2root, tree2array
import numpy as np
import ROOT as R



#%%

def npz_to_kintree(rfile):

    #rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/test.root")
    intree = rfile.Get('qfactortree')
    f = tree2array(intree)
    
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



    px_beam = f ['px_beam']
    py_beam = f['py_beam']
    pz_beam = f['pz_beam']
    e_beam = f['e_beam']

    
    
    px = []
    py = []
    pz = []
    e  = []
    
    
    e_bm = []
    px_bm = []
    py_bm = []
    pz_bm = []
    
    for i in range(len(px_pr)):
        #x components
        ux = px_pr[i]
        vx = px_etapr[i]
        wx = px_pi0[i]
        
        #y components
        uy = py_pr[i]
        vy = py_etapr[i]
        wy = py_pi0[i]
    
        #z components
        uz = pz_pr[i]
        vz = pz_etapr[i]
        wz = pz_pi0[i]
    
        #energy
        ue = e_pr[i]
        ve = e_etapr[i]
        we = e_pi0[i]
        
        #4 vector combonents for beam
        be = e_beam[i]
        bx = px_beam[i]
        by = py_beam[i]
        bz = pz_beam[i]
    
        #x, y, z  and e of final states combined 
        ax = np.array([ux, vx, wx], dtype = 'float64')
        ay = np.array([uy, vy, wy], dtype = 'float64')
        az = np.array([uz, vz, wz], dtype = 'float64')
        ae = np.array([ue, ve, we], dtype = 'float64')
        
        
        #x, y, z  and e of beam 
        be = np.array([be], dtype = 'float64')
        bx = np.array([bx], dtype = 'float64')
        by = np.array([by], dtype = 'float64')
        bz = np.array([bz], dtype = 'float64')
        
        
        px.append(ax)
        py.append(ay)
        pz.append(az)
        e.append(ae)
        
        e_bm.append(bz)
        px_bm.append(bx)
        py_bm.append(by)
        pz_bm.append(bz)
    return e, px, py, pz, e_bm, px_bm, py_bm,  pz_bm


#%%

rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/test.root")
npz_to_kintree(rfile)
e, px, py, pz, e_bm, px_bm, py_bm, pz_bm = npz_to_kintree(rfile)

#%%

e = np.array(e)
px = np.array(px)
py = np.array(py)
pz = np.array(pz)


e_bm = np.array(e_bm)
px_bm = np.array(px_bm)
py_bm = np.array(py_bm)
pz_bm = np.array(pz_bm)

#%%

e.dtype = [('E_FinalState', np.float64, (3,))]
px.dtype = [('Px_FinalState', np.float64, (3,))]
py.dtype = [('Py_FinalState', np.float64, (3,))]
pz.dtype = [('Pz_FinalState', np.float64, (3,))]


e_bm.dtype = [('E_Beam', np.float64, (1,))]
px_bm.dtype = [('Px_Beam', np.float64, (1,))]
py_bm.dtype = [('Py_Beam', np.float64, (1,))]
pz_bm.dtype = [('Pz_Beam', np.float64, (1,))]

#for test 
qs_weights = np.random.random_sample((px.shape[0]))
qs_weights.dtype = [('signal_weights', np.float64, (1,))]


#%%
array2root(e, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/PWA.root", treename='kin', mode = 'recreate')

array2root(px, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')

array2root(py, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(pz, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')


array2root(e_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(px_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(py_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(pz_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')

array2root(qs_weights, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')


#%%
'''
rf = R.TFile("for_PWA.root")
intree = rf.Get('kin')
d = tree2array(intree)

'''
    

