#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:50:09 2018

@author: ricktjwong
"""

import numpy as np
from scipy import signal
import numba
import time


# Inputs
a = np.random.rand(100,100)

# Convert to numpy array
arr = np.asarray(a,float)    

# Define kernel for convolution                                         
kernel = np.array([[0,1,0],
                   [1,0,1],
                   [0,1,0]]) 


@numba.njit
def convolve():
    out = np.zeros(a.shape)
    for r in range(1, a.shape[0] - 1):
        for c in range(1, a.shape[1] - 1):
            value = kernel * a[(r - 1):(r + 2), (c - 1):(c + 2)]
            out[r, c] = value.sum()/4
    return out


start = time.time()
# Perform 2D convolution with input data and kernel 
out = signal.convolve2d(arr, kernel, boundary='wrap', mode='same')/kernel.sum()
end = time.time()
print(end - start)
print(out[1:-1,1:-1].sum())


start = time.time()
out = convolve()
end = time.time()
print(end - start)
print(out[1:-1,1:-1].sum())
