#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:50:09 2018

@author: ricktjwong
"""

import numpy as np
import time

# Inputs
a = np.random.rand(200,200)

def jacobi_loop(a):
    out = np.zeros(a.shape)
    for r in range(1, a.shape[0] - 1):
        for c in range(1, a.shape[1] - 1):
            out[r, c] = (a[r][c-1] + a[r][c+1] + a[r-1][c] + a[r+1][c]) / 4
    return out

def jacobi_solver(x):
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

start = time.time()
out = jacobi_loop(a)
end = time.time()
print(end - start)
print(out[1:-1,1:-1].sum())

start = time.time()
out = jacobi_solver(a).copy()
end = time.time()
print(end - start)
print(out[1:-1,1:-1].sum())
