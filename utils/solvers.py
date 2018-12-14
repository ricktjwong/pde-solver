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
    
    y[self.m_idx_y1:self.m_idx_y2,
      self.m_idx_x1:self.m_idx_x2] += self.rho
           
    return y


def jacobi_loop(self, x):
    """
    Inefficient implementation of Jacobi method which is O(n^2), looping
    through the entire array
    """
    out = x.copy()
    for r in range(1, x.shape[0] - 1):
        for c in range(1, x.shape[1] - 1):
            out[r, c] = (x[r][c-1] + x[r][c+1] + x[r-1][c] + x[r+1][c]) / 4
    out[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] += self.rho
    return out


def jacobi_convolve(self, x):
    """
    Efficient Jacobi solver which uses scipy's convolution to operate on the
    entire array at once
    """
    # Define kernel for convolution                                         
    kernel = np.array([[0,1,0],
                       [1,0,1],
                       [0,1,0]])
    out = signal.convolve2d(x, kernel, boundary='wrap', mode='same') / \
    kernel.sum()
    out[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] += self.rho
    return out


def gauss_seidel(self, x):
    # Update matrix for the fins        
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
    for i in range(self.m_idx_y1, self.m_idx_y2):
        for j in range(self.m_idx_x1, self.m_idx_x2):        
            x[i][j] = 1/4 * (x[i-1][j] + x[i+1][j] 
                            + x[i][j-1] + x[i][j+1]) + self.rho
    return x


def SOR_solver(self, x):    
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
    # Update matrix for the microprocessor
    for i in range(self.m_idx_y1, self.m_idx_y2):
        for j in range(self.m_idx_x1, self.m_idx_x2):        
            x[i][j] = (1 - w) * x[i][j] + w/4 * \
                      (x[i-1][j] + x[i+1][j] + 
                       x[i][j-1] + x[i][j+1] + 4 * self.rho)
    return x

   
def SOR_solver2(self, x):
    w = 0.8
    backup = x.copy()
    x_u = np.roll(x, -1, 0)
    x_u[-1] = np.zeros((1, x.shape[1]))
    
    x_d = np.roll(x, 1, 0)
    x_d[0] = np.zeros((1, x.shape[1]))
    
    x_r = np.roll(x, 1, 1)
    x_r[:,0] = np.zeros((1, x.shape[0]))
    
    x_l = np.roll(x, -1, 1)
    x_l[:,-1] = np.zeros((1, x.shape[0]))

    b = np.ones(x.shape)
    r = np.ones(x.shape)
    b[::2,::2] = 0
    b[1::2,1::2] = 0
    r -= b
    
    y = (1 - w) * x + w * (x_d + x_u + x_r + x_l) / 4
    y[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] += w * self.rho
    final_red = b * y + r * backup
    
    backup = final_red.copy()
    x_u = np.roll(final_red, -1, 0)
    x_u[-1] = np.zeros((1, final_red.shape[1]))
    
    x_d = np.roll(final_red, 1, 0)
    x_d[0] = np.zeros((1, final_red.shape[1]))
    
    x_r = np.roll(final_red, 1, 1)
    x_r[:,0] = np.zeros((1, final_red.shape[0]))
    
    x_l = np.roll(final_red, -1, 1)
    x_l[:,-1] = np.zeros((1, final_red.shape[0]))
    
    y = (1 - w) * final_red + w * (x_d + x_u + x_r + x_l) / 4
    y[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] += w * self.rho
    final_black = r * y + b * backup
#    print (final_black)
    return final_black
     