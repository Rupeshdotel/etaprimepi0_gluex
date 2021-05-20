#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 18:40:24 2021

@author: rupeshdotel

select unique events based on pair (run number, event number) 
"""

import numpy as np



#%%
def remove_duplicates(array):
    
    na = []
    for a in array:
        if a not in na:
            na.append(a)
    return na


def find_duplicates(np_array):
    # np_array is a numpy array
    # return a list with np_array values that contain no duplicates and a
    # corresponding list of indices pointing to identical values in np_array
    # optionall return a list of elements of duplicates into array
    
    array_no_duplicates = np.array(remove_duplicates(np_array))
    # sort the array 
    #array_no_duplicates.sort()
    #index_list = List()
    index_list = []
    #value_index_list = []
    for a in array_no_duplicates:
        # compare arrays element wise
        np_array = np.array(np_array)
        comp = np_array == a
        # find where all elements of a are equal to the same set in np_array
        
        if len(np_array.shape) > 1:
            equal = np.apply_along_axis(np.min, 1, comp)
        else:
            equal = comp
        # get the indices
        in_dup = np.where(equal)[0]
        #in_dup = List()
        #index_list.append(List(in_dup)) 
        index_list.append( list(in_dup) )
        #value_index_list.append(list(in_dup)[0])
    #if get_value_index:
     #   return array_no_duplicates, index_list, value_index_list
    #else:
    return array_no_duplicates, index_list


    
#%%
#d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluex_17_18_cons_incamo.npz')
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_cons.npz')
en = d['event_num']
rn = d['run_num']
pol = d['pol']

#use m variables and cost for data only to do q factors
mpi0 = d['mpi0']
meta = d['meta']
metap = d['metap']
metappi0 = d['metappi0']
#cost_etap = d['cost_etap']

#cost_etap_gj = d['cost_etap_gj']
#phi_etap_gj = d['phi_etap_gj']

cost_etap_gj = d['cost_etap']
phi_etap_gj = d['phi_etap']
#mant = d['mant']


#2) required for both MC and data

px_pr = d['px_pr']
px_etapr = d['px_etapr']
px_pi0 = d['px_pi0']

py_pr = d['py_pr']
py_etapr = d['py_etapr']
py_pi0 = d['py_pi0']

pz_pr = d['pz_pr']
pz_etapr = d['pz_etapr']
pz_pi0 = d['pz_pi0']

e_pr = d['e_pr']
e_etapr = d['e_etapr']
e_pi0 = d['e_pi0']

px_beam = d['px_beam']
py_beam = d['py_beam']
pz_beam = d['pz_beam']
e_beam = d['e_beam']

#%%


#for test
#en = en[:500]
#rn = rn[:500]
#mpi0 = mpi0[:500]

pair = np.array([rn, en]) #run number event number pair
pair = pair.T


pair = pair.tolist() # get the pair as list outside the function 

#%%
pair_nd, il = find_duplicates(pair[:]) #pair_nd is pair with duplicates removed

#%%
len_il = [len(i) for i in il]
len_il = np.array(len_il)


sel = (len_il == 1) 
unique_pair = pair_nd[sel] #unique pair is where an event number shows up only once


unique_en = unique_pair[:,1] # select corresponding event numbers

#%%

# create a mask to select other variables  for unique events
mask = []
for i in en:
    if i in unique_en:
        mask.append(True)
    else:
        mask.append(False)

mask = np.array(mask)

#%%
#get different variables corresponding to unique events variables
mpi0 = mpi0[mask]
meta = meta[mask]
metap = metap[mask]
metappi0 = metappi0[mask]
cost_etap_gj = cost_etap_gj[mask]
phi_etap_gj = phi_etap_gj[mask]
#mant = mant[mask]


#2) required for both MC and data

px_pr = px_pr[mask]
px_etapr = px_etapr[mask]
px_pi0 = px_pi0[mask]

py_pr = py_pr[mask]
py_etapr = py_etapr[mask]
py_pi0 = py_pi0[mask]

pz_pr = pz_pr[mask]
pz_etapr = pz_etapr[mask]
pz_pi0 = pz_pi0[mask]

e_pr = e_pr[mask]
e_etapr = e_etapr[mask]
e_pi0 = e_pi0[mask]

px_beam = px_beam[mask]
py_beam = py_beam[mask]
pz_beam = pz_beam[mask]
e_beam = e_beam[mask]
pol = pol[mask]

#%%

np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/gluexI_unique_cons.npz',
         
            #1) for data only
           
            #use m variables and cost for data only to do q factors
            mpi0 = mpi0,
            meta = meta,
            metap = metap,
            metappi0 = metappi0,
            cost_etap_gj = cost_etap_gj,
            phi_etap_gj = phi_etap_gj,
            #mant = mant,
            
            #2) required for both MC and data
            
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
            pol = pol
           
            
            
           )


        
            
            





    