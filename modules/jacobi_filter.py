#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:50:09 2018

@author: ricktjwong
"""

import numpy as np
from scipy import signal

# Inputs
a = [[1,2,3,4],[9,8,7,6],[1,4,5,8],[7,8,9,10]]

# Convert to numpy array
arr = np.asarray(a,float)    

# Define kernel for convolution                                         
kernel = np.array([[0,1,0],
                   [1,0,1],
                   [0,1,0]]) 

# Perform 2D convolution with input data and kernel 
out = signal.convolve2d(arr, kernel, boundary='wrap', mode='same')/kernel.sum()

print(np.array(a))
print(np.array(out))
