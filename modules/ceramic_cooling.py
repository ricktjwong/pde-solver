#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 17:22:14 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

h = 1E-5                    # Step size h (in mm)
k_m = 0.15                  # Conductivity of silicon Microchip in W/mm K
k_c = 0.23                  # Conductivity of ceramic block in W/mm K
del_T = 0.5 / k_m           # Change in temperature of microprocessor every s
T_a = 20                    # Ambient temperature
alpha_m = h * 2.62 / k_m    # Constant for natural convection
alpha_c = h * 2.62 / k_c    # Constant for natural convection
rows = 5
cols = 16


def update_boundaries(mesh):
    rows = mesh.shape[0]
    cols = mesh.shape[1]
    c_l_mesh = mesh[1, 0:3]
    m_l_mesh = mesh[2:-1, 0:3]
    c_r_mesh = mesh[1, cols-3:]
    m_r_mesh = mesh[2:-1, cols-3:]
    c_u_mesh = mesh[0:3, 1:-1]
    m_d_mesh = mesh[rows-3:, 1:-1]
    
    c_l_mesh[0] = c_l_mesh[-1] - alpha_c * (c_l_mesh[1] - T_a) ** (4/3)
    m_l_mesh[:,0] = m_l_mesh[:,-1] - alpha_m * (m_l_mesh[:,1] - T_a) ** (4/3)
    c_r_mesh[-1] = c_r_mesh[0] - alpha_c * (c_r_mesh[1] - T_a) ** (4/3)
    m_r_mesh[:,-1] = m_r_mesh[:,0] - alpha_m * (m_r_mesh[:,1] - T_a) ** (4/3)
    c_u_mesh[0] = c_u_mesh[-1] - alpha_c * (c_u_mesh[1] - T_a) ** (4/3)
    m_d_mesh[-1] = m_d_mesh[0] - alpha_m * (m_d_mesh[1] - T_a) ** (4/3)
    
    mesh[1, 0:3] = c_l_mesh
    mesh[2:-1, 0:3] = m_l_mesh
    mesh[1, cols-3:] = c_r_mesh
    mesh[2:-1, cols-3:] = m_r_mesh
    mesh[0:3, 1:-1] = c_u_mesh
    mesh[rows-3:, 1:-1] = m_d_mesh
    return mesh

mesh = np.zeros((rows, cols))
mesh[1, 1:-1] = 20
mesh[2:-1, 1:-1] = 20      # Initialise the inner mesh with T > T_a
print(mesh)

mesh = update_boundaries(mesh).copy()
print(mesh)

all_mesh = []

n = 0
while (n < 20000):
    update = mesh.copy()
    for j in range(1, cols-1):
        for i in range(1, rows-1):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    update[2:-1, 1:-1] += del_T
#    if np.linalg.norm(update[1:-1, 1:-1]/mesh[1:-1, 1:-1] - 1) < 1E-4:
#        break
    update = update_boundaries(update).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1

print(mesh)

x = []
m = []
for i in all_mesh:
#    x.append(i[1][1])
    m.append(i[2][1])
y = [i for i in range(len(m))]
plt.figure(1)
#plt.plot(y, x, "--", c='r')
plt.plot(y, m, "--", c='r')

plt.figure(2)
mesh_min = mesh[1:-1, 1:-1].min() - 0.01
mesh_max = mesh[1:-1, 1:-1].max()
plt.imshow(mesh, vmin=mesh_min, vmax=mesh_max)
