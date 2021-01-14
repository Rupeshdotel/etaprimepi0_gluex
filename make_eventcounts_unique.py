#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:55:41 2020

@author: rupeshdotel
"""

import numpy as np

#%%

# necessary tools to remove duplicate entries
#------------------------------------------------------
def remove_duplicates(array):
    na = []
    for a in array:   # try : if not found  go to except, if found skip except and go to return
        try:
            na.index(a)
        except:
            na.append(a)
    return na

#----------------------------------------------------------------------
def find_duplicates(np_array, get_value_index = False):
    # np_array is a numpy array
    # return a list with np_array values that contain no duplicates and a
    # corresponding list of indices pointing to identical values in np_array
    # optionally return a list of elements of duplicates into array
    array_no_duplicates = remove_duplicates( list(np_array) )
    # sort the array 
    array_no_duplicates.sort()
    index_list = []
    value_index_list = []
    for a in array_no_duplicates:
        in_dup = np.where(np_array == a)[0]
        index_list.append(list(in_dup))
         #choice 1 : select the first combo 
        #value_index_list.append(list(in_dup)[0])
        # choice 2: select the second/last combo from the multicombo events
        if len(list(in_dup)) == 1:
            value_index_list.append(list(in_dup)[0])
        #else:
         #    value_index_list.append(list(in_dup)[-1])
        
        
    if get_value_index:
        return array_no_duplicates, index_list, value_index_list
    else:
        return array_no_duplicates, index_list
    
#%%

#f = np.load('/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_17/qfactortree/new_eventsprompt_18_re.npz')
#f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/new_eventsprompt_17.npz')

# read the selected events, the event selection is already done
f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/selected_events.npz')
#i = 100
event_num = f['event_num'][:]
kinfit_CL = f['kinfit_CL'][:]
chisq_ndf = f['chisq_ndf'][:]
num_combos = f['num_combos'][:]
combo_number = f['combo_number'][:]
mpi0 = f['mpi0'][:]
meta = f['meta'][:]
metap = f['metap'][:]

metappi0 = f['metappi0'][:]

mpi013 = f['mpi013'][:]
mpi024 = f['mpi024'][:]
mpi014 = f['mpi014'][:]
mpi023 = f['mpi023'][:]

mant = f['mant'][:]
num_unusedshowers = f['num_unusedshowers'][:]

mpipp = f['mpipp'][:]
mpi0p = f['mpi0p'][:]
mpippimpi0 = f['mpippimpi0'][:]

pi0costhetaGJ = f['cost_pi0'][:]
pi0phiGJ = f['pi0phiGJ'][:]

etaprimecosthetaGJ = f['cost_etap'][:]
etaprimephiGJ = f['etaprimephiGJ'][:]

photon1_sq = f['photon1_sq'][:]
photon2_sq = f['photon2_sq'][:]
photon3_sq = f['photon3_sq'][:]
photon4_sq = f['photon4_sq'][:]

#%%
en ,il, vil = find_duplicates(event_num, True)
en = np.array(en)
il = np.array(il, dtype = object)
vil = np.array(vil)
#%%

event_num = event_num[vil]
kinfit_CL = kinfit_CL[vil]
chisq_ndf = chisq_ndf[vil]
num_combos = num_combos[vil]
combo_number = combo_number[vil]
mpi0 = mpi0[vil]
meta = meta[vil]
metap = metap[vil]

metappi0 = metappi0[vil]

mpi013 = mpi013[vil]
mpi024 = mpi024[vil]
mpi014 = mpi014[vil]
mpi023 = mpi023[vil]

mant = mant[vil]
num_unusedshowers = num_unusedshowers[vil]

mpipp = mpipp[vil]
mpi0p = mpi0p[vil]
mpippimpi0 = mpippimpi0[vil]

pi0costhetaGJ = pi0costhetaGJ[vil]
pi0phiGJ = pi0phiGJ[vil]

etaprimecosthetaGJ = etaprimecosthetaGJ[vil]
etaprimephiGJ = etaprimephiGJ[vil]

photon1_sq = photon1_sq[vil]
photon2_sq = photon2_sq[vil]
photon3_sq = photon3_sq[vil]
photon4_sq = photon4_sq[vil]


#%%
#np.savez('/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_17/qfactortree/unique_new_eventspromt_18_re.npz',
         
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/unique_selected_events.npz',         
            event_num = event_num,
            kinfit_CL = kinfit_CL,
            chisq_ndf = chisq_ndf,
            num_combos = num_combos,
            combo_number = combo_number,
            mpi0 = mpi0,
            meta = meta,
            metap = metap,
            
            metappi0 = metappi0,
            
            mpi013 = mpi013,
            mpi024 = mpi024,
            mpi014 = mpi014,
            mpi023 = mpi023,
            
            mant = mant,
            num_unusedshowers = num_unusedshowers,
            
            mpipp = mpipp,
            mpi0p = mpi0p,
            mpippimpi0 = mpippimpi0,
            
            cost_pi0 = pi0costhetaGJ,
            pi0phiGJ = pi0phiGJ,
            
            cost_etap = etaprimecosthetaGJ,
            etaprimephiGJ = etaprimephiGJ,
            
            photon1_sq = photon1_sq,
            photon2_sq = photon2_sq,
            photon3_sq = photon3_sq,
            photon4_sq = photon4_sq
            
            )












