#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 22:58:42 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

h = 0.05                    # Step size h (in mm)
k_m = 0.150                 # Conductivity of silicon Microchip in W/mm K
k_c = 0.230                 # Conductivity of ceramic block in W/mm K
k_a = 0.248                 # Conductivity of aluminium fin in W/mm K
del_T = h ** 2 * 0.5 / k_m  # Change in temperature of microprocessor every s
T_a = 20                    # Ambient temperature
alpha_m = h * 2.62 / k_m    # Constant for natural convection for silicon
alpha_c = h * 2.62 / k_c    # Constant for natural convection for ceramic
alpha_a = h * 2.62 / k_a    # Constant for natural convection for aluminium
rows = 37
cols = 22


def setup():
    mesh = np.zeros((rows, cols))
    mesh[1:31, 1:3] = T_a           # Initialise fin 1 to T_a
    mesh[1:31, 7:9] = T_a           # Initialise fin 2 to T_a
    mesh[1:31, 13:15] = T_a         # Initialise fin 3 to T_a
    mesh[1:31, 19:21] = T_a         # Initialise fin 4 to T_a
    mesh[31:35, 1:-1] = T_a         # Initialise the ceramic block at T_a
    mesh[35][4:18] = T_a            # Initialise the microprocessor at T_a
    return mesh


def update_boundaries(m, c):
    rows = m.shape[0]
    cols = m.shape[1]
    l_mesh = m[1:-1, 0:3]
    r_mesh = m[1:-1, cols-3:]
    u_mesh = m[0:3, 1:-1]
    d_mesh = m[rows-3:, 1:-1]

#    all_mesh = [l_mesh, r_mesh, u_mesh, d_mesh]
    for i in range(len(l_mesh[:,1])):
        if l_mesh[:,1][i] < 20:
            l_mesh[:,1][i] = 20
        else:
            l_mesh[i][0] = l_mesh[i][-1] - c * (l_mesh[i][1] - T_a) ** (4/3)
    for i in range(len(r_mesh[:,1])):
        if r_mesh[:,1][i] < 20:
            r_mesh[:,1][i] = 20
        else:
            r_mesh[i][-1] = r_mesh[i][0] - c * (r_mesh[i][1] - T_a) ** (4/3)
    for i in range(len(u_mesh[1])):
        if u_mesh[1][i] < 20:
            u_mesh[1][i] = 20
        else:
            u_mesh[:,i][0] = u_mesh[:,i][-1] - c * (u_mesh[:,i][1] - T_a) ** (4/3)
    for i in range(len(d_mesh[1])):
        if d_mesh[1][i] < 20:
            d_mesh[1][i] = 20    
        else:
            d_mesh[:,i][-1] = d_mesh[:,i][0] - c * (d_mesh[:,i][1] - T_a) ** (4/3)
#    for k in all_mesh:
#        for j in range(len(k[0])):
#            for i in range(len(k)):
#                if k[i][j] < 20:
#                    print(k[i][j])
#                    k[i][j] = 20

#    l_mesh[:,0] = l_mesh[:,-1] - c * (l_mesh[:,1] - T_a) ** (4/3)
#    r_mesh[:,-1] = r_mesh[:,0] - c * (r_mesh[:,1] - T_a) ** (4/3)
#    u_mesh[0] = u_mesh[-1] - c * (u_mesh[1] - T_a) ** (4/3)
#    d_mesh[-1] = d_mesh[0] - c * (d_mesh[1] - T_a) ** (4/3)
    
    m[1:-1, 0:3] = l_mesh
    m[1:-1, cols-3:] = r_mesh
    m[0:3, 1:-1] = u_mesh
    m[rows-3:, 1:-1] = d_mesh
    

def update_all_components(m):
    c_mesh = m[30:36].copy()
    update_boundaries(c_mesh, alpha_c)
    m_mesh = m[33:, 3:19].copy()
    update_boundaries(m_mesh, alpha_m)
    # 4 fins of width 2px, height 30, spacing 4
    f_mesh = [m[0:32, 0 + 6*i: 4 + 6*i].copy() for i in range(4)]
    for f in f_mesh:
        update_boundaries(f, alpha_a)
    return c_mesh, m_mesh, f_mesh


def update_mesh(m, c_mesh, m_mesh, f_mesh):
    m[30:36] = c_mesh.copy()
    for i in range(4):
        m[0:31, 0 + 6*i: 4 + 6*i] = f_mesh[i][0:-1].copy()
    m[35:, 4:18] = m_mesh[2:, 1:-1].copy()
    return m

mesh = setup()
c_mesh, m_mesh, f_mesh = update_all_components(mesh)
mesh = update_mesh(mesh, c_mesh, m_mesh, f_mesh).copy()

all_mesh = []

n = 0
while (n < 1000):
    update = mesh.copy()
    # Update matrix for the ceramic block
    for j in range(1, cols-1):
        for i in range(31, 35):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    # Update matrix for the microprocessor
    r = 35
    for j in range(4, 18):
        update[r][j] = 1/4 * (mesh[r-1][j] + mesh[r+1][j] 
                            + mesh[r][j-1] + mesh[r][j+1])
    # Update matrix for the fins        
    for k in range(4):
        for j in range(1 + 6*k, 3 + 6*k):
            for i in range(1, 31):
                update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                    + mesh[i][j-1] + mesh[i][j+1])
        
    update[35][4:18] += del_T
#    if np.linalg.norm(update[1:-1, 1:-1]/mesh[1:-1, 1:-1] - 1) < 1E-4:
#        break
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
