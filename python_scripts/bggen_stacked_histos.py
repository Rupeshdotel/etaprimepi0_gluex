#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 15:39:53 2020

@author: rupeshdotel

read root histos and  and convert to python histos and plot in stack
"""

import ROOT as R
import bin_info2 as bi
import LT.box as B
import matplotlib.pyplot as plt



#%%
#read the root file
rfile  = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/bggen_etaprimepi0_KL_largerthan_0.001_mm2cut_coh_be_intime_etaprimewindow.root")


#%%
#convert root 1d histo to python 1d histo
def get_1d_root_histo(r_histo):
    # get approriate histogram
    hr = bi.get_histo_data(r_histo)
    # make arrays with the correct shape
    hr.make_arrays()
    # make a histo1d
    hd = B.histo(bin_center = hr.xb, bin_content = hr.cont, bin_error = hr.dcont ) 
    # transfer the labels and title
    title = r_histo.GetName()
    xl = (r_histo.GetXaxis()).GetTitle() 
    hd.xlabel = xl
    hd.title = title
    return hd
#%%

#histos for different topologies ( nomenclature comes from root histos   )
h_ep_pi0 = get_1d_root_histo(rfile.metappi0)
h_pi0_eta_etap = get_1d_root_histo(rfile.hInvariantMass_ThrownTopology_2)
h_pi0_eta_omega = get_1d_root_histo(rfile.hInvariantMass_ThrownTopology_3)
h_pi0_eta = get_1d_root_histo(rfile.hInvariantMass_ThrownTopology_0)
h_2pi0_eta = get_1d_root_histo(rfile.hInvariantMass_ThrownTopology_1)
h_2pi0_omega = get_1d_root_histo(rfile.hInvariantMass_ThrownTopology_4)


#%%
#draw the stakced histos
plt.figure();
h_ep_pi0.plot(filled  = False)
h_pi0_eta_etap.plot(filled = True, color = 'green', alpha = 0.6)
h_pi0_eta_omega.plot(filled  = True, color = 'black', alpha = 0.6)
h_pi0_eta.plot(filled  = True, color = 'blue', alpha = 0.6)
h_2pi0_eta.plot(filled  = True, color = "red", alpha = 0.6)
h_2pi0_omega.plot(filled  = True, color = "yellow", alpha = 0.6)
plt.title("thrown topology")
#plt.legend(['thrown topology'])
plt.legend([("Total reconstructed $\eta^{'}\pi^{0}$ invariant mass  "),("$4\gamma\pi^{+}\pi^{-}p[\pi^{0},\eta,\eta^{'}]$"), 
            ("$4\gamma\pi^{+}\pi^{-}p[\pi^{0},\eta,\omega]$"), 
            ("$4\gamma\pi^{+}\pi^{-}p[\pi^{0},\eta]$"),
           ("$4\gamma\pi^{+}\pi^{-}p[2\pi^{0},\eta]$"),
           ("$4\gamma\pi^{+}\pi^{-}p[2\pi^{0},\omega]$")])










