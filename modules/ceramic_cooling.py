#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 03:09:37 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time

scale = 2
h = 1E-3 / scale            # Step size h (in m)
k_m = 150                   # Conductivity of silicon Microchip in W/m K
k_c = 230                   # Conductivity of ceramic block in W/m K
k_a = 248                   # Conductivity of aluminium fin in W/m K
rho = 0.25 * h ** 2 * 500 * 1E6 / k_m  # Change in temperature of microprocessor every s
T_a = 20                    # Ambient temperature
alpha_m = h * 2.62 / k_m    # Constant for natural convection for silicon
alpha_c = h * 2.62 / k_c    # Constant for natural convection for ceramic
rows = 3 * scale + 2
cols = 20 * scale + 2

# Initialise row and column indices for ceramic block
c_idx_x1, c_idx_x2 = 1, cols - 1
c_idx_y1, c_idx_y2 = 1, 2*scale + 1
# Initialise row and column indices for microprocessor block
m_idx_x1 = int((cols - 14*scale)/2)
m_idx_x2 = m_idx_x1 + 14 * scale
m_idx_y1, m_idx_y2 = c_idx_y2, c_idx_y2 + 1 * scale


def setup():
    mesh = np.zeros((rows, cols))        
    mesh[c_idx_y1:c_idx_y2, c_idx_x1:c_idx_x2] = T_a                  # Initialise the ceramic block at T_a
    mesh[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] = T_a                  # Initialise the microprocessor at T_a
    return mesh


def update_boundaries(m, c):
    rows = m.shape[0]
    cols = m.shape[1]
    l_mesh = m[1:-1, 0:3]
    r_mesh = m[1:-1, cols-3:]
    u_mesh = m[0:3, 1:-1]
    d_mesh = m[rows-3:, 1:-1]

    l_mesh[:,0] = l_mesh[:,-1] - c * (l_mesh[:,1] - T_a) ** (4/3)
    r_mesh[:,-1] = r_mesh[:,0] - c * (r_mesh[:,1] - T_a) ** (4/3)
    u_mesh[0] = u_mesh[-1] - c * (u_mesh[1] - T_a) ** (4/3)
    d_mesh[-1] = d_mesh[0] - c * (d_mesh[1] - T_a) ** (4/3)
    
    m[1:-1, 0:3] = l_mesh
    m[1:-1, cols-3:] = r_mesh
    m[0:3, 1:-1] = u_mesh
    m[rows-3:, 1:-1] = d_mesh
    

def update_all_components(m):
    c_mesh = m[c_idx_y1-1:c_idx_y2+1, c_idx_x1-1:c_idx_x2+1].copy()
    update_boundaries(c_mesh, alpha_c)
    m_mesh = m[m_idx_y1-2:m_idx_y2+1, m_idx_x1-1:m_idx_x2+1].copy()
    update_boundaries(m_mesh, alpha_m)
    return c_mesh, m_mesh


def update_mesh(m, c_mesh, m_mesh):
    m[c_idx_y1-1:c_idx_y2+1, c_idx_x1-1:c_idx_x2+1] = c_mesh.copy()
    m[m_idx_y1:m_idx_y2+1, m_idx_x1-1:m_idx_x2+1] = m_mesh[2:].copy()
    return m

mesh = setup()
plt.imshow(mesh)
c_mesh, m_mesh = update_all_components(mesh)
mesh = update_mesh(mesh, c_mesh, m_mesh).copy()

all_mesh = []
plt.imshow(mesh)

start = time.time()
n = 0
while (True):
    update = mesh.copy()
    # Define kernel for convolution                                         
    kernel = np.array([[0,1,0],
                       [1,0,1],
                       [0,1,0]]) 
    # Perform 2D convolution with input data and kernel 
    out = signal.convolve2d(mesh, kernel,
                            boundary='wrap', mode='same')/kernel.sum()
    # Update matrix for the ceramic block    
    update[c_idx_y1:c_idx_y2, c_idx_x1:c_idx_x2] = \
    out[c_idx_y1:c_idx_y2, c_idx_x1:c_idx_x2]
    # Update matrix for the microprocessor
    update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] = \
    out[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]
    update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] += rho
    if np.mean(update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]) / \
       np.mean(mesh[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]) - 1 < 1E-6:
        break
    c_mesh, m_mesh = update_all_components(update)
    update = update_mesh(update, c_mesh, m_mesh).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1
end = time.time()
print(end - start)

x = []
for i in all_mesh:
    x.append(np.mean(i[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]))
y = [i for i in range(len(x))]

plt.figure(2)
plt.plot(y, x, "--", c='r')

plt.figure(3)
masked = np.ma.masked_where(mesh < 0.01, mesh)
plt.imshow(masked, cmap="rainbow")
plt.show()
