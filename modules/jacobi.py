#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 22:03:45 2018

@author: ricktjwong
"""

import numpy as np

def jacobi_method(A, b):
    """
    Solve the matrix equation Ax = b using the iterative Jacobi method
    Ax = b => (D + R)x = b
    Dx + Rx = b
    """
    rows = A.shape[0]
    cols = A.shape[1]

    # Initialise empty arrays for the diagonal and inverse diagonal matrices
    D = np.array([[0 for i in range(cols)]
                 for i in range(rows)], dtype=float)

    inv_D = np.array([[0 for i in range(cols)]
                     for i in range(rows)], dtype=float)

    # Update diagonal and inverse diagonal matrices
    for j in range(cols):
        for i in range(rows):
            if i == j:
                D[i][i] = A[i][i]
                inv_D[i][i] = 1 / A[i][i]

    R = A - D
    # Initialise the iteration with x_0 = [1, 1, 1]
    x_0 = np.array([1, 1, 1], dtype=float)

    # Iterate until the change in new x value is smaller than 1E-14
    while True:
        x_1 = np.dot(inv_D, (b - np.dot(R, x_0)))
        if np.linalg.norm(x_1 / x_0 - 1) < 1E-14:
            break
        x_0 = x_1

    return x_0
