#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:32:17 2021

@author: rupeshdotel

generate kin tree for PWA from clean samples in the form of numpy arrays

"""

from root_numpy import array2tree, array2root, tree2array
import numpy as np
#import ROOT as R



#%%

def npz_to_kintree(filename):

    
    #intree = filename.Get('qfactortree')
    #f = tree2array(intree)
    f = filename
    
    
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
        ax = np.array([ux, vx, wx], dtype = 'float32')
        ay = np.array([uy, vy, wy], dtype = 'float32')
        az = np.array([uz, vz, wz], dtype = 'float32')
        ae = np.array([ue, ve, we], dtype = 'float32')
        
        
        #x, y, z  and e of beam 
        be = np.array([be], dtype = 'float32')
        bx = np.array([bx], dtype = 'float32')
        by = np.array([by], dtype = 'float32')
        bz = np.array([bz], dtype = 'float32')
        
        
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

#file = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortee/recon_mc.root")
#intree = file.Get('qfactortree')
#d = tree2array(intree)

file = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_pol_amo.npz')
#file = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events_recon_mc_pwa.npz')

#file = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI.npz')
e, px, py, pz, e_bm, px_bm, py_bm, pz_bm = npz_to_kintree(file)

#%%
#convert lists to arrays
e = np.array(e)
px = np.array(px)
py = np.array(py)
pz = np.array(pz)


e_bm = np.array(e_bm)
px_bm = np.array(px_bm)
py_bm = np.array(py_bm)
pz_bm = np.array(pz_bm)




#%%
#specify type explicitly, not doing so will create error while writing the tree
"""
e = e.astype([('E_FinalState', np.float32, (3,))])
px = px.astype([('Px_FinalState', np.float32, (3,))])
py = py.astype([('Py_FinalState', np.float32, (3,))])
pz = pz.astype([('Pz_FinalState', np.float32, (3,))])


e_bm = e_bm.astype([('E_Beam', np.float32, (1,))])
px_bm = px_bm.astype([('Px_Beam', np.float32, (1,))])
py_bm = py_bm.astype([('Py_Beam', np.float32, (1,))])
pz_bm = pz_bm.astype([('Pz_Beam', np.float32, (1,))])
"""

e.dtype = [('E_FinalState', np.float32, (3,))]
px.dtype = [('Px_FinalState', np.float32, (3,))]
py.dtype = [('Py_FinalState', np.float32, (3,))]
pz.dtype = [('Pz_FinalState', np.float32, (3,))]


e_bm.dtype = [('E_Beam', np.float32, (1,))]
px_bm.dtype = [('Px_Beam', np.float32, (1,))]
py_bm.dtype = [('Py_Beam', np.float32, (1,))]
pz_bm.dtype = [('Pz_Beam', np.float32, (1,))]


#for test 
#qs_weights = np.random.random_sample((px.shape[0]))


#%%
#num_fs = np.zeros_like(e_bm)
num_fs = np.repeat(3, e_bm.shape[0])
num_fs = num_fs.astype([('NumFinalState', np.int32, (1,))])


pol = file['pol']
pol.dtype = [('Polarization_Angle', np.int32, (1,))]
#num_fs.dtype = [('NumFinalState', np.int, (1,))]


#%%
"""
#only for data
#read qvalues from different file
f_w = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/weights.npz')
qs_weights = f_w['weights']

#convert qvalues to a float (originally its a double)
qs_weights = qs_weights.astype([('Weight', np.float32, (1,))])
#qs_weights.dtype = [('Weight', np.float32, (1,))]
"""
qs_weights = file['qf']
qs_weights = qs_weights.astype([('Weight', np.float32, (1,))])



#%%
#write the kin tree here, first one recreate, then update
array2root(num_fs, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'recreate')

array2root(e, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')

array2root(px, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')

array2root(py, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(pz, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')


array2root(e_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(px_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(py_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')
array2root(pz_bm, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')

array2root(pol, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')



array2root(qs_weights, filename = "/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/for_PWA.root", treename='kin', mode = 'update')


#%%
'''
rf = R.TFile("for_PWA.root")
intree = rf.Get('kin')
d = tree2array(intree)

'''
    

