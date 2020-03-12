# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 18:54:31 2018

@author: ifsw_linux2

utilities for handling stl data

"""

import numpy as np
import matplotlib.pyplot as plt

#### importieren einer ASCII STL file als punktewolke im numpy array format

def read_ASCII_STL(filename):
    """
    reads ASCII STL, returns np.array[n,3] of vertexes
    """
    f = open(filename)
    lines = f.readlines()
    lines = [l for l in lines if 'vertex' in l]
    data = np.asanyarray([[float (x) for x in [e for e in l.rstrip('\n').split(' ') if e != ''][1:]] for l in lines])
    f.close()
    return data
    
def delete_multiples(points):
    return np.unique(points, axis = 0)


#filename = '/home/ifsw_linux2/OpenFOAM/OF_main_v0_l/Bin/Scripts/pyGeom/tests/Surfaces/kapillare_neu.stl'
#read_ASCII_STL(filename)
#points = delete_multiples(read_ASCII_STL(filename))