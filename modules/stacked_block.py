#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:25:09 2018

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
cols = 22 


def update_boundaries(mesh):
    rows = mesh.shape[0]
    cols = mesh.shape[1]

    c_l_mesh = mesh[1:rows-2, 0:3]
    c_r_mesh = mesh[1:rows-2, cols-3:]
    c_u_mesh = mesh[0:3, 1:-1]    
    c_dl_mesh = mesh[rows-4:-1, 1:4]    # Slice to update bottom left of ceramic
    c_dr_mesh = mesh[rows-4:-1, 18:-1]  # Slice to update bottom right of ceramic
    m_dd_mesh = mesh[rows-3:, 4:18]     # Slice to update bottom of microprocessor
    
    c_l_mesh[:,0] = c_l_mesh[:,-1] - alpha_c * (c_l_mesh[:,1] - T_a) ** (4/3)
    c_r_mesh[:,-1] = c_r_mesh[:,0] - alpha_c * (c_r_mesh[:,1] - T_a) ** (4/3)
    c_u_mesh[0] = c_u_mesh[-1] - alpha_c * (c_u_mesh[1] - T_a) ** (4/3)
    c_dl_mesh[-1] = c_dl_mesh[0] - alpha_c * (c_dl_mesh[1] - T_a) ** (4/3)
    c_dr_mesh[-1] = c_dr_mesh[0] - alpha_c * (c_dr_mesh[1] - T_a) ** (4/3)        
    m_dd_mesh[-1] = m_dd_mesh[0] - alpha_m * (m_dd_mesh[1] - T_a) ** (4/3)
    
    mesh[1:rows-2, 0:3] = c_l_mesh
    mesh[1:rows-2, cols-3:] = c_r_mesh
    mesh[0:3, 1:-1] = c_u_mesh
    mesh[rows-4:-1, 1:4] = c_dl_mesh
    mesh[rows-4:-1, 18:-1] = c_dr_mesh
    mesh[rows-3:, 4:18] = m_dd_mesh
    return mesh

# Initialise the microprocessor and ceramic block to a start temperature of T_a
mesh = np.zeros((rows, cols))
mesh[3][4:18] = T_a
mesh[1:3][:, 1:-1] = T_a
                         
print(mesh)
mesh = update_boundaries(mesh).copy()
print(mesh)

all_mesh = []

n = 0
while (True):
    update = mesh.copy()
    # Update matrix for the ceramic block
    for j in range(1, cols-1):
        for i in range(1, rows-2):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    # Update matrix for the microprocessor
    for j in range(4, 18):
        update[3][j] = 1/4 * (mesh[2][j] + mesh[4][j] 
                            + mesh[3][j-1] + mesh[3][j+1])
        
    update[3][4:18] += del_T

    if np.linalg.norm(update[1:-1, 1:-1]/mesh[1:-1, 1:-1] - 1) < 1E-4:
        break
    update = update_boundaries(update).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1

print(mesh)

x = []
m = []
n = []
for i in all_mesh:
    x.append(i[1][1])
    m.append(i[1][2])
    n.append(i[3][4])
y = [i for i in range(len(x))]

plt.figure(1)
#plt.plot(y, x, "--", c='r')
#plt.plot(y, m, "--", c='b')
plt.plot(y, n, "--", c='g')

plt.figure(2)
mesh_min = mesh[1:-1, 1:-1].min() - 0.01
mesh_max = mesh[1:-1, 1:-1].max()
plt.imshow(mesh, vmin=mesh_min, vmax=mesh_max)
