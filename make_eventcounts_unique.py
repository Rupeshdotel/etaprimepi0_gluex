#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 19:33:18 2020

@author: rupeshdotel
"""

import ROOT as R
from root_numpy import tree2array
import matplotlib.pyplot as plt
import numpy as np
import LT.box as B

#%%
#get the input root file with tree
rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/qfactortree/qfactortree_0303.root")
#rfile = R.TFile("/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_mc/qfactortree/qfactortree_sim_prompt_flat.root")
#rfile = R.TFile("/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_17/qfactortree/qfactortree_17.root")
intree = rfile.Get('qfactortree')

#%%
array = tree2array(intree)


# get the variables for q_factor analysis
event_num = array['event_num']
kinfit_CL = array['kinfit_CL']
chisq_ndf = array['chisq_ndf']
num_combos = array['num_combos']
combo_number = array['combo_number']

# invariant masses
pi0mass = array['pi0mass']
etaprimemass = array['etaprimemass']
etaprimepi0mass = array['etaprimepi0mass']

#angular variables
pi0costhetaGJ = array['pi0costhetaGJ']
pi0phiGJ = array['pi0phiGJ']

etaprimecosthetaGJ = array['etaprimecosthetaGJ']
etaprimephiGJ = array['etaprimephiGJ']

pippimpi0 = array['pippimpi0']
pipp = array['pipp']
pi0p = array['pi0p']
dt = array['dt']
etamass = array['etamass']
BE = array['BEa']

num_combos = array['num_combos']
combo_number_unique = array['combo_number']

#%%
np.savez('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/eventsaccid_gluexI.npz', 
         #'/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_17/qfactortree/unique_eventsprompt_17'
         pi0mass = pi0mass , 
         etaprimemass= etaprimemass,
         etaprimepi0mass = etaprimepi0mass, 
         
         pi0costhetaGJ = pi0costhetaGJ,
         pi0phiGJ = pi0phiGJ, 
         etaprimecosthetaGJ = etaprimecosthetaGJ,
         etaprimephiGJ = etaprimephiGJ,
         
         event_num = event_num,
         kinfit_CL = kinfit_CL,
         chisq_ndf = chisq_ndf,
         
         pippimpi0 = pippimpi0,
         pipp = pipp,
         pi0p = pi0p,
         dt = dt,
         etamass = etamass,
         BE = BE,
         
         num_combos = num_combos,
         combo_number = combo_number
         
         
         )


#%%
'''plt.plot(event_num,kinfit_CL,".")
plt.plot(event_num,kinfit_CL)
plt.hist(num_combos, bins=20)
plt.hist(chisq_ndf, bins=20)
plot(chisq_ndf,kinfit_CL,".")

B.plot_exp(chisq_ndf,kinfit_CL, logy = True)
B.plot_exp(chisq_ndf,kinfit_CL)'''

#%%
#f = np.load('/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_18/qfactortree/eventsprompt_18_08.npz')
'''f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/eventsprompt_18_08.npz')

#i = 1006093
pi0mass = f['pi0mass'][:] 
etaprimemass= f['etaprimemass'][:]
etaprimepi0mass = f['etaprimepi0mass'][:] 

pi0costhetaGJ = f['pi0costhetaGJ'][:]
pi0phiGJ = f['pi0phiGJ'][:] 
etaprimecosthetaGJ = f['etaprimecosthetaGJ'][:]
etaprimephiGJ = f['etaprimephiGJ'][:]

num_combos = f['num_combos'][:]
combo_number = f['combo_number'][:]

event_num = f['event_num'][:]
kinfit_CL = f['kinfit_CL'][:]
chisq_ndf = f['chisq_ndf'][:]

pippimpi0 = f['pippimpi0'][:]
pipp = f['pipp'][:]
pi0p = f['pi0p'][:]
dt = f['dt'][:]
etamass = f['etamass'][:]
BE = f['BE'][:]'''



#%%
'''f = np.load('/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_18/qfactortree/eventsprompt_18.npz')
#f = np.load('/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/eventsprompt_17.npz')
pi0mass = f['pi0mass']
etaprimemass= f['etaprimemass']
etaprimepi0mass = f['etaprimepi0mass']

pi0costhetaGJ = f['pi0costhetaGJ']
pi0phiGJ = f['pi0phiGJ']
etaprimecosthetaGJ = f['etaprimecosthetaGJ']
etaprimephiGJ = f['etaprimephiGJ']

num_combos = f['num_combos']
combo_number = f['combo_number']

event_num = f['event_num']
kinfit_CL = f['kinfit_CL']
chisq_ndf = f['chisq_ndf']

pippimpi0 = f['pippimpi0']
pipp = f['pipp']
pi0p = f['pi0p']
dt = f['dt']
etamass = f['etamass']
BE = f['BE']
'''

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
#----------------------------------------------------------------------
#%%

#find_duplicates(event_num)
event_num_unique ,indexlist,value_index_list = find_duplicates(event_num, True)
event_num_unique = np.array(event_num_unique)
indexlist = np.array(indexlist)
valueindexlist = np.array(value_index_list)


pi0mass_unique = pi0mass[valueindexlist]
etaprimemass_unique = etaprimemass[valueindexlist]
etaprimepi0mass_unique = etaprimepi0mass[valueindexlist]

event_num_unique = event_num[valueindexlist]
kinfit_CL_unique = kinfit_CL[valueindexlist]
chisq_ndf_unique = chisq_ndf[valueindexlist]
num_combos_unique = num_combos[valueindexlist]
combo_number_unique = combo_number[valueindexlist]


pi0costhetaGJ_unique = pi0costhetaGJ[valueindexlist]
pi0phiGJ_unique = pi0phiGJ[valueindexlist]

etaprimecosthetaGJ_unique = etaprimecosthetaGJ[valueindexlist]
etaprimephiGJ_unique = etaprimephiGJ[value_index_list]

pippimpi0_unique = pippimpi0[valueindexlist]
pipp_unique = pipp[valueindexlist]
pi0p_unique = pi0p[valueindexlist]
dt_unique = dt[valueindexlist]
etamass_unique = etamass[valueindexlist]
BE_unique = BE[valueindexlist]

#print(pi0mass_unique)
#plt.hist(pi0mass_unique, bins = 120)

#%%
'''np.savez(#'/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_mc/unique_eventsprompt_17.npz', 
         '/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_17/qfactortree/unique_eventsprompt_17.npz'
         pi0mass_unique_sc = pi0mass_unique , 
         etaprimemass_unique_sc = etaprimemass_unique,
         etaprimepi0mass_unique_sc = etaprimepi0mass_unique, 
         
         pi0costhetaGJ_unique_sc = pi0costhetaGJ_unique,
         pi0phiGJ_unique_sc = pi0phiGJ_unique, 
         etaprimecosthetaGJ_unique_sc = etaprimecosthetaGJ_unique,
         etaprimephiGJ_unique_sc = etaprimephiGJ_unique,
         
         event_num_unique_sc = event_num_unique,
         kinfit_CL_unique_sc = kinfit_CL_unique,
         chisq_ndf_unique_sc = chisq_ndf_unique,
         
         pippimpi0_unique_sc = pippimpi0_unique,
         pipp_unique_sc = pipp_unique,
         pi0p_unique_sc = pi0p_unique,
         dt_unique_sc = dt_unique,
         etamass_unique_sc = etamass_unique,
         BE_sc = BE_unique
         
         )

'''

#%%
np.savez(#'/w/halld-scifs17exp/halld2/home/rupesh/halld/pi0etapr__B4_M35_M7_M17_18/qfactortree/unique_eventsprompt_18.npz',
         '/Users/rupeshdotel/analysis/work/pi0pippimeta/data/qfactor_data/uniqueventsaccid_0303.npz',
         pi0mass_unique  = pi0mass_unique, 
         etaprimemass_unique = etaprimemass_unique,
         etaprimepi0mass_unique = etaprimepi0mass_unique, 
         
         pi0costhetaGJ_unique = pi0costhetaGJ_unique,
         pi0phiGJ_unique = pi0phiGJ_unique, 
         etaprimecosthetaGJ_unique = etaprimecosthetaGJ_unique,
         etaprimephiGJ_unique = etaprimephiGJ_unique,
         
         event_num_unique = event_num_unique,
         kinfit_CL_unique = kinfit_CL_unique,
         chisq_ndf_unique = chisq_ndf_unique,
         
         pippimpi0_unique = pippimpi0_unique,
         pipp_unique = pipp_unique,
         pi0p_unique = pi0p_unique,
         dt_unique = dt_unique,
         etamass_unique = etamass_unique,
         BE_unique = BE_unique,
         num_combos_unique = num_combos_unique,
         combo_number_unique = combo_number_unique
         
         
         )



#%%
#------------weighted combo scheme for test ----------------
'''
def weighted_combo_events(np_array):
    weighted_combo_list = []
    single_combo_list = []
    nd = [len(x) for x in indexlist]
    nd = np.array(nd)   
    single_entry = nd == 1
    multiple_entry = nd > 1

    for i in indexlist[single_entry]:
        # only one combo entry per event  so no need to multiply by 1
        single_combo_list.append(list(np_array[i]))
        
#print("singlecombolist", single_combo_list)

    for i in indexlist[multiple_entry]:
        weight = 1/len(i) #weight by number of combos in the particular event
    #pi0mass[i]
        #print(np_array[i])
        weighted_combo_elements = [] 
        for j in np_array[i]:
       
            #print("ic", j)
            #print("wc", weight * j)
            wce = weight * j
            weighted_combo_elements.append(wce)
        #print("wcl", wcl)
        #print("added", sum(wcl))
        weighted_combo_list.append(sum(weighted_combo_elements))
#print("weightedcombolist", weighted_combo_list)

    flattened_singlecombo_list = [y for x in single_combo_list for y in x]
    weighted_events = flattened_singlecombo_list + weighted_combo_list
    #plt.hist(weighted_events, bins=30)
    return  weighted_events

#%%

weighted_combo_events(pi0mass) 


  

'''




    
'''

 for j in i:
        pi0mass[j] = weight*pi0mass[j]
    sum(pi0mass[j])
    weighted_combo.append(pi0mass[j])
print(weighted_combo)
en
len(en)
il
vil
len(vil)

plt.hist(pi0mass[vil])
combo_number[vil] #selects first combo of the event and rejects the  rest
en
il
nd = [len(x) for x in il]
nd
nd = np.array(nd)
nd
sel = nd > 1
en[sel]
il[sel]
il = np.array(il)
il[sel]
il[sel][0]
en[sel][0]
chisq_ndf[il[sel][0]]
len(en[sel])
chisq_ndf[il[sel][41]]
nd.max
nd.max()
h = B.histo(nd, range=(-0.5,10.5),bins = 11)
h.plot_exp()
num_combos[il[0]]
num_combos[il[1]]
num_combos[il[2]]
num_combos[il[2]].max
num_combos[il[2]].max()
num_combos.max()
event_num[num_combos==100]
chisq_ndf[num_combos==100]
kinfit_CL[num_combos==100]


#%%
#count the frequency : number of repetitions
a = np.array(event_num)
unique, counts = numpy.unique(a, return_counts=True)
unique, counts = np.unique(a, return_counts=True)
dict(zip(unique, counts))
'''






    