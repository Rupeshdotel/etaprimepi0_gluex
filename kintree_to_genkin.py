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
from ROOT import TVector3


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
    cos_theta = np.zeros_like(Px)
    
    #declare  TLorentzVector
    Beam_P4 = TLorentzVector()
    Proton_P4 = TLorentzVector()
    Etaprime_P4 = TLorentzVector()
    Pi0_P4 = TLorentzVector()
    Etaprimepi0_P4 = TLorentzVector()
    
    #declare stuff to convert to GJ frame
    boostGJ = TVector3()
    Beam_P4_GJ = TLorentzVector()
    Proton_P4_GJ = TLorentzVector()
    Etaprime_P4_GJ = TLorentzVector()
    Pi0_P4_GJ = TLorentzVector()
    Etaprimepi0_P4_GJ = TLorentzVector()
    
    z_GJ = TVector3()
    z_hat_GJ = TVector3()
    
    y_GJ = TVector3()
    y_hat_GJ = TVector3()
    
    x_hat_GJ = TVector3()
    vetaprime = TVector3() 
    
    #for i in range(len(Px)):
    for i in range(100):
        #m = np.sqrt(E[i][0]**2 - (Px[i][0]**2 + Py[i][0]**2 + Pz[i][0]**2))
        #metap[i] = m
        
        
        #get  TLorentzVector vectors
        Beam_P4.SetPxPyPzE(Px_beam[i], Py_beam[i], Pz_beam[i], E_beam[i])
        Proton_P4.SetPxPyPzE(Px[i][0], Py[i][0], Pz[i][0], E[i][0])
        Etaprime_P4.SetPxPyPzE(Px[i][1], Py[i][1], Pz[i][1], E[i][1])
        Pi0_P4.SetPxPyPzE(Px[i][2], Py[i][2], Pz[i][2], E[i][2])
        
        #do stuff here
        Etaprimepi0_P4 = Etaprime_P4 + Pi0_P4
        
        #get the boost vector
        boostGJ = (-1) * (Etaprimepi0_P4.Vect())*(1.0/Etaprimepi0_P4.E())
        
        
        Beam_P4_GJ = Beam_P4
        Proton_P4_GJ = Proton_P4
        Etaprime_P4_GJ = Etaprime_P4
        Pi0_P4_GJ = Pi0_P4
        Etaprimepi0_P4_GJ = Etaprimepi0_P4
        
        
        
        #boost vectors to GJ frame
        Beam_P4_GJ.Boost(boostGJ)
        Proton_P4_GJ.Boost(boostGJ)
        Etaprime_P4_GJ.Boost(boostGJ)
        Pi0_P4_GJ.Boost(boostGJ)
        Etaprimepi0_P4_GJ.Boost(boostGJ)
        
        z_GJ.SetXYZ(Beam_P4_GJ.X(),Beam_P4_GJ.Y(),Beam_P4_GJ.Z())
        z_hat_GJ = z_GJ.Unit() 
        
        y_GJ = Beam_P4.Vect().Cross(Etaprimepi0_P4.Vect())
        y_hat_GJ = y_GJ.Unit() 
        
        x_hat_GJ = y_hat_GJ.Cross(z_hat_GJ)
        
        vetaprime.SetXYZ(Etaprime_P4_GJ.Vect()*x_hat_GJ, Etaprime_P4_GJ.Vect()*y_hat_GJ, Etaprime_P4_GJ.Vect()*z_hat_GJ)
        cos_theta[i]  = vetaprime.CosTheta()
        
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
                    cos_theta = cos_theta,
                    beam_energy = beam_energy)

"""load the file as 
g = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/MC/gen_events.npz',  allow_pickle = True)
get numpy array as 
metaprimepi0 = g['metaprimepi0']"""
    
    
    
    
    
    