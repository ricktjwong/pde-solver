#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:17:59 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time

scale = 1
h = 1E-3 / scale            # Step size h (in m)
k_m = 150                   # Conductivity of silicon Microchip in W/m K
k_c = 230                   # Conductivity of ceramic block in W/m K
k_a = 248                   # Conductivity of aluminium fin in W/m K
rho = 0.25 * h ** 2 * 500 * 1E6 / k_m  # Change in temperature of microprocessor every s
T_a = 20                    # Ambient temperature
alpha_m = h * 2.62 / k_m    # Constant for natural convection for silicon
alpha_c = h * 2.62 / k_c    # Constant for natural convection for ceramic
alpha_a = h * 2.62 / k_a    # Constant for natural convection for aluminium
n_fins = 4
b = 4
c = 2
f_h = 30
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
    mesh = np.zeros((rows, cols))        
    for i in range(n_fins):
        mesh[f_idx_y1:f_idx_y2, (T*i*scale+1):(c+T*i)*scale+1] = T_a  # Initialise fin 1 to T_a
    mesh[fb_idx_y1:fb_idx_y2, fb_idx_x1:fb_idx_x2] = T_a              # Initialise the fin block to T_a        
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
    

def update_all_boundaries(m):
    c_mesh = m[c_idx_y1-1:c_idx_y2+1, c_idx_x1-1:c_idx_x2+1].copy()
    update_boundaries(c_mesh, alpha_c)
    m_mesh = m[m_idx_y1-2:m_idx_y2+1, m_idx_x1-1:m_idx_x2+1].copy()
    update_boundaries(m_mesh, alpha_m)
    fb_mesh = m[fb_idx_y1-1:fb_idx_y2+1, fb_idx_x1-1:fb_idx_x2+1].copy()
    update_boundaries(fb_mesh, alpha_a)
    # 4 fins of width 2px, height 30, spacing 4
    f_mesh = [m[f_idx_y1-1:f_idx_y2+1, (T*i*scale):(c+T*i)*scale+2].copy() \
                for i in range(n_fins)]
    for f in f_mesh:
        update_boundaries(f, alpha_a)
    return c_mesh, m_mesh, fb_mesh, f_mesh


def update_mesh(m, c_mesh, m_mesh, fb_mesh, f_mesh):
    m[fb_idx_y1-1:fb_idx_y2+1, fb_idx_x1-1:fb_idx_x2+1] = fb_mesh.copy()
    for i in range(n_fins):
        m[f_idx_y1-1:f_idx_y2, (T*i*scale):(c+T*i)*scale+2] = \
        f_mesh[i][0:-1].copy()
    m[c_idx_y1:c_idx_y2+1, c_idx_x1-1:c_idx_x2+1] = c_mesh[1:].copy()
    m[m_idx_y1:m_idx_y2+1, m_idx_x1-1:m_idx_x2+1] = m_mesh[2:].copy()
    return m

mesh = setup()
plt.imshow(mesh)
c_mesh, m_mesh, fb_mesh, f_mesh = update_all_boundaries(mesh)
mesh = update_mesh(mesh, c_mesh, m_mesh, fb_mesh, f_mesh).copy()

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
    update[fb_idx_y1:fb_idx_y2, fb_idx_x1:fb_idx_x2] = \
    out[fb_idx_y1:fb_idx_y2, fb_idx_x1:fb_idx_x2]
    # Update matrix for the fins
    for i in range(n_fins):
        update[f_idx_y1:f_idx_y2, (T*i*scale+1):(c+T*i)*scale+1] = \
        out[f_idx_y1:f_idx_y2, (T*i*scale+1):(c+T*i)*scale+1]
    update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2] += rho
    if np.mean(update[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]) / \
       np.mean(mesh[m_idx_y1:m_idx_y2, m_idx_x1:m_idx_x2]) - 1 < 1E-6:
        break
    c_mesh, m_mesh, fb_mesh, f_mesh = update_all_boundaries(update)
    update = update_mesh(update, c_mesh, m_mesh, fb_mesh, f_mesh).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1
end = time.time()
print(end - start)

x = []
for i in all_mesh:
    x.append(i[31][1])    
y = [i for i in range(len(x))]

plt.figure(2)
plt.plot(y, x, "--", c='r')

plt.figure(3)
masked = np.ma.masked_where(mesh < 0.01, mesh)
plt.imshow(masked, cmap="rainbow")
plt.show()
