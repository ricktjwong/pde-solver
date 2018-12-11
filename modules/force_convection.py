#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:55:58 2018
@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time

scale = 4
n_fins = 5
b = 5
c = 1
f_h = 54

h = 1E-3 / scale            # Step size h (in m)
k_m = 150                   # Conductivity of silicon Microchip in W/m K
k_c = 230                   # Conductivity of ceramic block in W/m K
k_a = 248                   # Conductivity of aluminium fin in W/m K
# Change in temperature of microprocessor every s
rho = 0.25 * h ** 2 * 500 * 1E6 / k_m
T_a = 20                    # Ambient temperature
# Constant for force convection for silicon
alpha_m = h * 2 / k_m * (11.54 + 5.7 * 20)
# Constant for force convection for ceramic
alpha_c = h * 2 / k_c * (11.54 + 5.7 * 20)
# Constant for force convection for aluminium
alpha_a = h * 2 / k_a * (11.54 + 5.7 * 20)
T = b+c
rows = (f_h + 7) * scale + 2
cols = (T*(n_fins-1) + c) * scale + 2

# Initialise row indices for fins
f_idx_y1, f_idx_y2 = 1, (f_h*scale)+1
# Initialise row and column indices for fin block
fb_idx_x1, fb_idx_x2 = 1, cols - 1
fb_idx_y1, fb_idx_y2 = f_idx_y2, f_idx_y2 + 4 * scale
# Initialise row and column indices for ceramic block
c_idx_x1 = int((cols - 20*scale)/2)
c_idx_x2 = c_idx_x1 + 20*scale
c_idx_y1, c_idx_y2 = fb_idx_y2, fb_idx_y2 + 2 * scale
# Initialise row and column indices for microprocessor block
m_idx_x1 = c_idx_x1 + 3 * scale
m_idx_x2 = m_idx_x1 + 14 * scale
m_idx_y1, m_idx_y2 = c_idx_y2, c_idx_y2 + 1 * scale


def setup():
    """
    Initialise the initial mesh with boundary points
    """
    mesh = np.zeros((rows, cols))
    # Initialise fin mesh points to T_a
    for i in range(n_fins):
        mesh[f_idx_y1:f_idx_y2, (T*i*scale+1):(c+T*i)*scale+1] = T_a
    # Initialise the fin block mesh points to T_a        
    mesh[fb_idx_y1:fb_idx_y2, fb_idx_x1:fb_idx_x2] = T_a
    # Initialise the ceramic block mesh points to T_a
    mesh[c_idx_y1:c_idx_y2, c_idx_x1:c_idx_x2] = T_a
    # Initialise the microprocessor mesh points to T_a
    mesh[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] = T_a
    return mesh


def update_boundaries(m, c):
    """
    Calculate the values for boundaries using the Central Difference Scheme
    (CDS). Each mesh (m) represents a rectangular component of the entire
    structure
    """
    rows = m.shape[0]
    cols = m.shape[1]
    # Extract the non-boundary mesh points for four sides of the block
    l_mesh = m[1:-1, 0:3]
    r_mesh = m[1:-1, cols-3:]
    u_mesh = m[0:3, 1:-1]
    d_mesh = m[rows-3:, 1:-1]
    # Update the boundary points using the CDS for four sides of block
    l_mesh[:,0] = l_mesh[:,-1] - c * (l_mesh[:,1] - T_a)
    r_mesh[:,-1] = r_mesh[:,0] - c * (r_mesh[:,1] - T_a)
    u_mesh[0] = u_mesh[-1] - c * (u_mesh[1] - T_a)
    d_mesh[-1] = d_mesh[0] - c * (d_mesh[1] - T_a)
    # Update the original mesh with the new boundary values
    m[1:-1, 0:3] = l_mesh
    m[1:-1, cols-3:] = r_mesh
    m[0:3, 1:-1] = u_mesh
    m[rows-3:, 1:-1] = d_mesh
    

def update_all_boundaries(m):
    """
    Update boundaries for each building block of the entire structure
    """
    c_mesh = m[c_idx_y1-1:c_idx_y2+1, c_idx_x1-1:c_idx_x2+1].copy()
    update_boundaries(c_mesh, alpha_c)
    m_mesh = m[m_idx_y1-2:m_idx_y2+1, m_idx_x1-1:m_idx_x2+1].copy()
    update_boundaries(m_mesh, alpha_m)
    fb_mesh = m[fb_idx_y1-1:fb_idx_y2+1, fb_idx_x1-1:fb_idx_x2+1].copy()
    update_boundaries(fb_mesh, alpha_a)
    f_mesh = [m[f_idx_y1-1:f_idx_y2+1, (T*i*scale):(c+T*i)*scale+2].copy() \
                for i in range(n_fins)]
    for f in f_mesh:
        update_boundaries(f, alpha_a)
    return c_mesh, m_mesh, fb_mesh, f_mesh


def update_mesh(m, c_mesh, m_mesh, fb_mesh, f_mesh):
    """
    Update original mesh with the boundary values of different components. The
    components are stacked in an order that preserves the correct boundaries.
    In this case we stack in the order: fin block, fins, ceramic, microchip
    """
    # Fin block with all boundaries
    m[fb_idx_y1-1:fb_idx_y2+1, fb_idx_x1-1:fb_idx_x2+1] = fb_mesh.copy()
    # Fins without the bottom boundary
    for i in range(n_fins):
        m[f_idx_y1-1:f_idx_y2, (T*i*scale):(c+T*i)*scale+2] = \
        f_mesh[i][0:-1].copy()
    # Ceramic without the top boundary
    m[c_idx_y1:c_idx_y2+1, c_idx_x1-1:c_idx_x2+1] = c_mesh[1:].copy()
    # Microprocessor without the top boundary
    m[m_idx_y1:m_idx_y2+1, m_idx_x1-1:m_idx_x2+1] = m_mesh[2:].copy()
    return m


def update_nonboundaries(m, n):
    """
    Update all the mesh points for the structure excluding the boundary values.
    m is the original mesh, n is the new mesh to be returned
    """
    # Update matrix for the ceramic block    
    n[c_idx_y1:c_idx_y2, c_idx_x1:c_idx_x2] = \
    m[c_idx_y1:c_idx_y2, c_idx_x1:c_idx_x2]
    # Update matrix for the microprocessor
    n[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] = \
    m[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]
    # Update matrix for the fin block
    n[fb_idx_y1:fb_idx_y2, fb_idx_x1:fb_idx_x2] = \
    m[fb_idx_y1:fb_idx_y2, fb_idx_x1:fb_idx_x2]
    # Update matrix for the fins
    for i in range(n_fins):
        n[f_idx_y1:f_idx_y2, (T*i*scale+1):(c+T*i)*scale+1] = \
        m[f_idx_y1:f_idx_y2, (T*i*scale+1):(c+T*i)*scale+1]


def solve_mesh(mesh, conv_ratio):
    """
    Iterate and solve for the temperature at every mesh point till the change
    in average temperature of the microprocessor is below the convergence_ratio
    """
    n = 0
    all_mesh = []
    while (True):
        update = mesh.copy()
        # Define Jacobi filter kernel for convolution                                         
        kernel = np.array([[0,1,0],
                           [1,0,1],
                           [0,1,0]]) 
        # Perform 2D convolution with input data and Jacobi filter kernel
        out = signal.convolve2d(mesh, kernel,
                                boundary='wrap', mode='same')/kernel.sum()
        update_nonboundaries(out, update)
        update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] += rho
        m_mean_temp_1 = np.mean(mesh[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2])
        m_mean_temp_2 = np.mean(update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2])
        if m_mean_temp_2 / m_mean_temp_1 - 1 < conv_ratio: break
        c_mesh, m_mesh, fb_mesh, f_mesh = update_all_boundaries(update)
        update = update_mesh(update, c_mesh, m_mesh, fb_mesh, f_mesh).copy()
        mesh = update.copy()
        all_mesh.append(mesh)
        n += 1
    return all_mesh


#mesh = setup()
#c_mesh, m_mesh, fb_mesh, f_mesh = update_all_boundaries(mesh)
#mesh = update_mesh(mesh, c_mesh, m_mesh, fb_mesh, f_mesh).copy()
#
#start = time.time()
#all_mesh = solve_mesh(mesh, 1E-6)
#end = time.time()
#print(end - start)
#
#x = []
#for i in all_mesh:
#    x.append(np.mean(i[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]))
#y = [i for i in range(len(x))]
#
#plt.figure(2)
#plt.plot(y, x, "--", c='r')
#
#plt.figure(3)
#final_mesh = np.zeros((rows, cols))
#update_nonboundaries(all_mesh[-1], final_mesh)
#
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#plt.imshow(masked, cmap="rainbow")
#plt.show()
