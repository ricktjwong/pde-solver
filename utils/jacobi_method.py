#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 16:26:11 2018

@author: ricktjwong
"""

import numpy as np

def jacobi_solver(x):
    """
    Takes an array x, applies the Jacobi pictorial operator and returns the
    new matrix y
    """
    x_u = np.roll(x, -1, 0)
    x_u[-1] = np.zeros((1, x.shape[1]))
    
    x_d = np.roll(x, 1, 0)
    x_d[0] = np.zeros((1, x.shape[1]))
    
    x_r = np.roll(x, 1, 1)
    x_r[:,0] = np.zeros((1, x.shape[0]))
    
    x_l = np.roll(x, -1, 1)
    x_l[:,-1] = np.zeros((1, x.shape[0]))
    
    y = (x_d + x_u + x_r + x_l) / 4
    return y
