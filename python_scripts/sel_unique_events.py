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
d = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/mm.npz')
en = d['en']
rn = d['rn']
mpi0 = d['mpi0']

#for test
#en = en[:500]
#rn = rn[:500]
#mpi0 = mpi0[:500]


pair = np.array([rn, en]) #run number event number pair
pair = pair.T


pair = pair.tolist() # get the pair as list outside the function 

pair_nd, il = find_duplicates(pair[:]) #pair_nd is pair with duplicates removed

len_il = [len(i) for i in il]
len_il = np.array(len_il)


sel = (len_il == 1) 
unique_pair = pair_nd[sel] #unique pair is where an event number shows up only once


unique_en = unique_pair[:,1] # select corresponding event numbers


# create a mask to select other variables  for unique events
mask = []
for i in en:
    if i in unique_en:
        mask.append(True)
    else:
        mask.append(False)

mask = np.array(mask)

#get different variables corresponding to unique events variables
mpi0 = mpi0[mask]







        
            
            





    