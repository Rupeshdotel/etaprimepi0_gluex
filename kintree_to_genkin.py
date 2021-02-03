#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 09:32:43 2021

@author: rupeshdotel

this script takes the kin tree as input (usually for generated MC ) and saves the desired 
kinematics into gen_events.npz file

particle order in final states 0:proton, 1:etaprime, 2:pi0

"""

import numpy as np
import ROOT as R
from root_numpy import tree2array
from ROOT import TLorentzVector

#%%
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/genamp_gen_resonance.root")

#%%
def get_gen_kin(rootfile):
    
    intree = rfile.Get('kin')
    d = tree2array(intree)
    
    Px = d['Px_FinalState']
    Py = d['Py_FinalState']
    Pz = d['Pz_FinalState']
    E = d['E_FinalState']

    Px_beam = d['Px_Beam']
    Py_beam = d['Py_Beam']
    Pz_beam = d['Pz_Beam']
    E_beam = d['E_Beam']
    
    #get numpy arrays to save later
    mproton = np.zeros_like(Px)
    metap = np.zeros_like(Px)
    mpi0 = np.zeros_like(Px)
    metaprimepi0 = np.zeros_like(Px)
    beam_energy = np.zeros_like(Px)
    
    #declare  TLorentzVector
    Beam_P4 = TLorentzVector()
    Proton_P4 = TLorentzVector()
    Etaprime_P4 = TLorentzVector()
    Pi0_P4 = TLorentzVector()
    Etaprimepi0_P4 = TLorentzVector()
    
    
    
    for i in range(len(Px)):
    #for i in range(100):
        #m = np.sqrt(E[i][0]**2 - (Px[i][0]**2 + Py[i][0]**2 + Pz[i][0]**2))
        #metap[i] = m
        
        
        #get  TLorentzVector vectors
        Beam_P4.SetPxPyPzE(Px_beam[i], Py_beam[i], Pz_beam[i], E_beam[i])
        Proton_P4.SetPxPyPzE(Px[i][0], Py[i][0], Pz[i][0], E[i][0])
        Etaprime_P4.SetPxPyPzE(Px[i][1], Py[i][1], Pz[i][1], E[i][1])
        Pi0_P4.SetPxPyPzE(Px[i][2], Py[i][2], Pz[i][2], E[i][2])
        
        
        Etaprimepi0_P4 = Etaprime_P4 + Pi0_P4
        
        
        mproton[i] = Proton_P4.M()
        metap[i] = Etaprime_P4.M()
        mpi0[i] = Pi0_P4.M()
        metaprimepi0[i] = Etaprimepi0_P4.M()
        beam_energy[i] = Beam_P4.E()
        
    
    return np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/gen_events.npz', 
                    mproton = mproton,
                    metap = metap, 
                    mpi0 = mpi0, 
                    metaprimepi0 = metaprimepi0, 
                    beam_energy = beam_energy)
    #return mproton, metap, mpi0, metaprimepi0, beam_energy
    
    
    
    
    
    