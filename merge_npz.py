#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 18:32:04 2021

@author: rupeshdotel
"""

import numpy as np

#%%

"""  

merge S17, S18, and Fall18 

example to concatenate different arrays of different files to a single array and store in a new file
a = np.array([1,2])
b = np.array([3,4])

np.savez('blah.npz', a = a)
np.savez('bla.npz', b = b)

f1 = np.load('blah.npz')
f2 = np.load('bla.npz')

c = np.concatenate([f1['a'], f2['b']])
np.savez('blahbla.npz', c = c )

"""