#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 03:40:36 2018

@author: ricktjwong
"""

import time
import modules.heat_structure as hst
import modules.utils.solvers as solv


def test_jacobi():
    """
    Assert that the three methods of Jacobi implementation return the same
    inner mesh values, and that the vectorised implementation is always faster
    that the loop version
    """
    hs = hst.HeatStructure(1, 2, 2, 5, 6, convection_type="forced",
                           solver=solv.jacobi_loop, conv_ratio=1E-3)
    start = time.time()
    jacobi_loop_T, n = hs.solve_mesh()
    end = time.time()
    jacobi_loop_time = (end - start)

    hs = hst.HeatStructure(1, 2, 2, 5, 6, convection_type="forced",
                           solver=solv.jacobi_solver, conv_ratio=1E-3)
    start = time.time()
    jacobi_fast_T, n = hs.solve_mesh()
    end = time.time()
    jacobi_fast_time = (end - start)
    
    hs = hst.HeatStructure(1, 2, 2, 5, 6, convection_type="forced",
                       solver=solv.jacobi_convolve, conv_ratio=1E-3)

    jacobi_conv_T, n = hs.solve_mesh()

    assert(jacobi_fast_time < jacobi_loop_time)
    assert(abs(jacobi_loop_T - jacobi_fast_T) < 1E-6)
    assert(abs(jacobi_loop_T - jacobi_conv_T) < 1E-6)
    assert(abs(jacobi_conv_T - jacobi_fast_T) < 1E-6)
