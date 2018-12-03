#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 00:04:17 2018

@author: ricktjwong
"""

import numpy as np
import modules.jacobi as jb

A = np.array([[1, 0, 0],
              [2, 3, 0],
              [0, 0, 5]], dtype=float)

x = np.array([1, 2, 3], dtype=float)
b = np.dot(A, x)

def test_jacobi_solver():
    x_0 = jb.jacobi_method(A, b)
    np.testing.assert_array_equal(x, x_0)
