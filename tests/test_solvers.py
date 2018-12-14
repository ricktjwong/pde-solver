#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 03:40:36 2018

@author: ricktjwong
"""

import numpy as np
import time
import utils.solvers as solvers

np.random.seed(1)
x = np.random.rand(500, 500)

def test_jacobi():
    """
    Assert that the three methods of Jacobi implementation return the same
    inner mesh values, and that the vectorised implementation is always faster
    that the loop version
    """
    start = time.time()
    out = solvers.jacobi_loop(x)
    end = time.time()
    jacobi_loop_time = (end - start)
    jacobi_loop = out[1:-1,1:-1].sum()

    start = time.time()
    out = solvers.jacobi_solver(x)
    end = time.time()
    jacobi_fast_time = (end - start)
    jacobi_fast = out[1:-1,1:-1].sum()
    
    out = solvers.jacobi_convolve(x)
    jacobi_conv = out[1:-1,1:-1].sum()    

    assert(jacobi_fast_time < jacobi_loop_time)
    assert(jacobi_loop == jacobi_fast == jacobi_conv)
