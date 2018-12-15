#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 16:26:11 2018

@author: ricktjwong
"""

import numpy as np
from scipy import signal


def jacobi_solver(self, x):
    """
    Efficient implementation of the Jacobi method. Takes an array x, applies
    the Jacobi pictorial operator and returns the new matrix y. This is done
    by array shifting using np.roll, improving the time complexity
    """
    # Shift array up, add zeros to the bottom row
    x_u = np.roll(x, -1, 0)
    x_u[-1] = np.zeros((1, x.shape[1]))
    # Shift array down, add zeros to the top row    
    x_d = np.roll(x, 1, 0)
    x_d[0] = np.zeros((1, x.shape[1]))
    # Shift array right, add zeros to the left column
    x_r = np.roll(x, 1, 1)
    x_r[:,0] = np.zeros((1, x.shape[0]))
    # Shift array left, add zeros to the right column
    x_l = np.roll(x, -1, 1)
    x_l[:,-1] = np.zeros((1, x.shape[0]))
    y = (x_d + x_u + x_r + x_l) / 4
    # Add source term to the microprocessor mesh. We cannot decouble the solver
    # logic from the source update logic, since the rest requires source term
    # updates within each iteration
    y[self.m_idx_y1:self.m_idx_y2,
      self.m_idx_x1:self.m_idx_x2] += self.rho
    return y


def jacobi_convolve(self, x):
    """
    Efficient Jacobi solver which uses scipy's convolution to operate on the
    entire array at once
    """
    # Define kernel for convolution - this is basically the Jacobi pictorial
    # operator                          
    kernel = np.array([[0,1,0],
                       [1,0,1],
                       [0,1,0]])
    # Perform a two dimensional convolution using scipy's convolve2d method.
    # This works because the kernel defines which terms contribute to the
    # resulting convolution. Averaging the result simply gets the Jacobi filter
    out = signal.convolve2d(x, kernel,
                            boundary='wrap', mode='same') / kernel.sum()
    # Add source term to the microprocessor mesh
    out[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] += self.rho
    return out


def jacobi_loop(self, x):
    """
    Inefficient implementation of Jacobi method which is O(n^2), looping
    through the entire array
    """
    out = x.copy()
    # Inefficient O(n^2) looping through of all cells of the mesh and updating
    # according to the Jacobi operator
    for r in range(1, x.shape[0] - 1):
        for c in range(1, x.shape[1] - 1):
            out[r, c] = (x[r][c-1] + x[r][c+1] + x[r-1][c] + x[r+1][c]) / 4
    out[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] += self.rho
    return out


def gauss_seidel(self, x):
    """
    Gauss Seidel iterative method, which is basically the Jacobi method, but
    the new values calculated after each update are used to calculate the next
    update, unlike the Jacobi method which generates a new array
    """
    # Update matrix for the fins
    # This loops through the the columns in each row
    for k in range(self.n_fins):
        for i in range(self.f_idx_y1, self.f_idx_y2):
            for j in range(self.T*k*self.scale+1,
                           (self.c+self.T*k)*self.scale+1):
                x[i][j] = 1/4 * (x[i-1][j] + x[i+1][j]
                                 + x[i][j-1] + x[i][j+1])
    # Update matrix for the fin block
    for i in range(self.fb_idx_y1, self.fb_idx_y2):
        for j in range(self.fb_idx_x1, self.fb_idx_x2):
            x[i][j] = 1/4 * (x[i-1][j] + x[i+1][j]
                            + x[i][j-1] + x[i][j+1])
    # Update matrix for the ceramic block
    for i in range(self.c_idx_y1, self.c_idx_y2):
        for j in range(self.c_idx_x1, self.c_idx_x2):
            x[i][j] = 1/4 * (x[i-1][j] + x[i+1][j]
                            + x[i][j-1] + x[i][j+1])
    # Update matrix for the microprocessor
    # Add source term to each microprocessor cell during each iteration
    for i in range(self.m_idx_y1, self.m_idx_y2):
        for j in range(self.m_idx_x1, self.m_idx_x2):
            x[i][j] = 1/4 * (x[i-1][j] + x[i+1][j]
                            + x[i][j-1] + x[i][j+1]) + self.rho
    return x


def SOR_solver(self, x):
    """
    Successive over-relaxation (SOR) method. This can be thought of as a
    combination between Gauss-Seidel and Jacobi method
    """
    w = 1.7
    # Update matrix for the fins
    for k in range(self.n_fins):
        for i in range(self.f_idx_y1, self.f_idx_y2):
            for j in range(self.T*k*self.scale+1,
                           (self.c+self.T*k)*self.scale+1):
                x[i][j] = (1 - w) * x[i][j] + w/4 * \
                          (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1])
    # Update matrix for the fin block
    for i in range(self.fb_idx_y1, self.fb_idx_y2):
        for j in range(self.fb_idx_x1, self.fb_idx_x2):
            x[i][j] = (1 - w) * x[i][j] + w/4 * \
                      (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1])
    # Update matrix for the ceramic block
    for i in range(self.c_idx_y1, self.c_idx_y2):
        for j in range(self.c_idx_x1, self.c_idx_x2):
            x[i][j] = (1 - w) * x[i][j] + w/4 * \
                      (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1])
    # Update matrix for the microprocessor. Add source term to each
    # microprocessor cell during each iteration, weighted by w
    for i in range(self.m_idx_y1, self.m_idx_y2):
        for j in range(self.m_idx_x1, self.m_idx_x2): 
            x[i][j] = (1 - w) * x[i][j] + w/4 * \
                      (x[i-1][j] + x[i+1][j] +
                       x[i][j-1] + x[i][j+1] + 4 * self.rho)
    return x
