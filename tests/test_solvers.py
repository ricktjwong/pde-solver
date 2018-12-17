#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 03:40:36 2018

@author: ricktjwong
"""

import sys
sys.path.append("../")
import numpy as np
import time
import modules.heat_structure as hst
import modules.utils.solvers as solv


def test_jacobi_shift():
    """
    Test the Jacobi array shifting method gives the average value of the
    neighbouring numbers
    """
    x = np.array([[1,2,3],
                  [4,5,6],
                  [7,10,12]])
    y = solv.jacobi_shift(x)
    assert(y[1][1] == 5.5)


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
    assert(abs(jacobi_loop_T[-1] - jacobi_fast_T[-1]) < 1E-6)
    assert(abs(jacobi_loop_T[-1] - jacobi_conv_T[-1]) < 1E-6)
    assert(abs(jacobi_conv_T[-1] - jacobi_fast_T[-1]) < 1E-6)
    
    
def test_all_solvers():
    """
    Assert that the three methods of Jacobi implementation return the same
    inner mesh values, and that the vectorised implementation is always faster
    that the loop version
    """
    hs = hst.HeatStructure(1, 2, 2, 5, 6, convection_type="forced",
                           solver=solv.SOR_solver, conv_ratio=1E-3)
    hs.solve_mesh()

    hs = hst.HeatStructure(1, 2, 2, 5, 6, convection_type="forced",
                           solver=solv.gauss_seidel, conv_ratio=1E-3)
    hs.solve_mesh()

    hs = hst.HeatStructure(1, 2, 2, 5, 6, convection_type="forced",
                       solver=solv.red_black_SOR, conv_ratio=1E-3)
    hs.solve_mesh()
    