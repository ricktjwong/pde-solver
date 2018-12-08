#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 01:53:55 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

h = 1E-3                                    # Step size h (in m)
k_m = 150                                   # Conductivity of silicon Microchip in W/m K
k_c = 230                                   # Conductivity of ceramic block in W/m K
k_a = 248                                   # Conductivity of aluminium fin in W/m K
rho = 0.25 * h ** 2 * 500 * 1E6 / k_m       # Change in temperature of microprocessor every s
T_a = 20                                    # Ambient temperature
alpha_m = h * 2 / k_m * (11.54 + 5.7 * 20)  # Constant for force convection for silicon
alpha_c = h * 2 / k_c * (11.54 + 5.7 * 20)  # Constant for force convection for ceramic
alpha_a = h * 2 / k_a * (11.54 + 5.7 * 20)  # Constant for force convection for aluminium
scale = 1
rows = 35 * scale + 2
cols = 20 * scale + 2


def setup():
    mesh = np.zeros((rows, cols))
    mesh[1:(30*scale)+1, 1:(2*scale)+1] = T_a                     # Initialise fin 1 to T_a
    mesh[1:(30*scale)+1, (6*scale)+1:(8*scale)+1] = T_a           # Initialise fin 2 to T_a
    mesh[1:(30*scale)+1, (12*scale)+1:(14*scale)+1] = T_a         # Initialise fin 3 to T_a
    mesh[1:(30*scale)+1, (18*scale)+1:(20*scale)+1] = T_a         # Initialise fin 4 to T_a
    mesh[(30*scale)+1:(34*scale)+1, 1:-1] = T_a                   # Initialise the ceramic block at T_a
    mesh[34*scale+1:35*scale+1, 3*scale+1:17*scale+1] = T_a       # Initialise the microprocessor at T_a
    return mesh


def update_boundaries(m, c):
    rows = m.shape[0]
    cols = m.shape[1]
    l_mesh = m[1:-1, 0:3]
    r_mesh = m[1:-1, cols-3:]
    u_mesh = m[0:3, 1:-1]
    d_mesh = m[rows-3:, 1:-1]

    l_mesh[:,0] = l_mesh[:,-1] - c * (l_mesh[:,1] - T_a)
    r_mesh[:,-1] = r_mesh[:,0] - c * (r_mesh[:,1] - T_a)
    u_mesh[0] = u_mesh[-1] - c * (u_mesh[1] - T_a)
    d_mesh[-1] = d_mesh[0] - c * (d_mesh[1] - T_a)
    
    m[1:-1, 0:3] = l_mesh
    m[1:-1, cols-3:] = r_mesh
    m[0:3, 1:-1] = u_mesh
    m[rows-3:, 1:-1] = d_mesh
    

def update_all_components(m):
    c_mesh = m[30*scale:34*scale+2].copy()
    update_boundaries(c_mesh, alpha_c)
    m_mesh = m[34*scale-1:, 3*scale:17*scale+2].copy()
    update_boundaries(m_mesh, alpha_m)
    # 4 fins of width 2px, height 30, spacing 4
    f_mesh = [m[0:30*scale+2, 6*i*scale:(2+6*i)*scale+2].copy() for i in range(4)]
    for f in f_mesh:
        update_boundaries(f, alpha_a)
    return c_mesh, m_mesh, f_mesh


def update_mesh(m, c_mesh, m_mesh, f_mesh):
    m[30*scale:34*scale+2] = c_mesh.copy()
    for i in range(4):
        m[0:30*scale+1, 6*i*scale:(2+6*i)*scale+2] = f_mesh[i][0:-1].copy()
    m[34*scale+1:, 3*scale:17*scale+2] = m_mesh[2:].copy()
    return m

mesh = setup()
plt.imshow(mesh)
c_mesh, m_mesh, f_mesh = update_all_components(mesh)
mesh = update_mesh(mesh, c_mesh, m_mesh, f_mesh).copy()

all_mesh = []
plt.imshow(mesh)

n = 0
while (True):
    update = mesh.copy()
    # Update matrix for the ceramic block
    for j in range(1, cols-1):
        for i in range(30*scale+1, 34*scale+1):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    # Update matrix for the microprocessor
    for j in range(3*scale+1, 17*scale+1):
        for i in range(34*scale+1, 35*scale+1):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    # Update matrix for the fins        
    for k in range(4):
        for j in range(1 + 6*k*scale, (2+6*k)*scale + 1):
            for i in range(1, 30*scale+1):
                update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                    + mesh[i][j-1] + mesh[i][j+1])
    update[34*scale+1:35*scale+1, 3*scale+1:17*scale+1] += rho
    if np.mean(update[34*scale+1:35*scale+1, 3*scale+1:17*scale+1]) / \
       np.mean(mesh[34*scale+1:35*scale+1, 3*scale+1:17*scale+1]) - 1 < 1E-6:
        break
    c_mesh, m_mesh, f_mesh = update_all_components(update)
    update = update_mesh(update, c_mesh, m_mesh, f_mesh).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1

x = []
for i in all_mesh:
    x.append(i[31][1])    
y = [i for i in range(len(x))]

plt.figure(2)
plt.plot(y, x, "--", c='r')

plt.figure(3)
masked = np.ma.masked_where(mesh < 0.01, mesh)
plt.imshow(masked, cmap="rainbow")
