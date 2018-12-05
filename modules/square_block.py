#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 00:47:27 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

h = 1E-5                    # Step size h (in mm)
k = 0.150                   # Conductivity of silicon Microchip in W/mm K
T_a = 20                    # Ambient temperature
alpha = h * 2.62 / k        # Constant for natural convection
rows = 10
cols = 10

def update_boundaries(mesh):
    rows = mesh.shape[0]
    cols = mesh.shape[1]
    l_mesh = mesh[1:-1, 0:3]
    r_mesh = mesh[1:-1, cols-3:]
    u_mesh = mesh[0:3, 1:-1]
    d_mesh = mesh[rows-3:, 1:-1]
    
    l_mesh[:,0] = l_mesh[:,-1] - alpha * (l_mesh[:,1] - T_a) ** (4/3)
    r_mesh[:,-1] = r_mesh[:,0] - alpha * (r_mesh[:,1] - T_a) ** (4/3)
    u_mesh[0] = u_mesh[-1] - alpha * (u_mesh[1] - T_a) ** (4/3)
    d_mesh[-1] = d_mesh[0] - alpha * (d_mesh[1] - T_a) ** (4/3)
    
    mesh[1:-1, 0:3] = l_mesh
    mesh[1:-1, cols-3:] = r_mesh
    mesh[0:3, 1:-1] = u_mesh
    mesh[rows-3:, 1:-1] = d_mesh
    return mesh

mesh = np.zeros((rows, cols))
mesh[1:-1, 1:-1] = 100      # Initialise the inner mesh with T > T_a
mesh = update_boundaries(mesh).copy()

print(mesh)

all_mesh = []

n = 0
while (True):
    update = mesh.copy()
    for j in range(1, cols-1):
        for i in range(1, rows-1):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    if np.linalg.norm(update[1:-1, 1:-1]/mesh[1:-1, 1:-1] - 1) < 1E-4:
        break
    update = update_boundaries(update).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1

print(mesh)

x = []
m = []
for i in all_mesh:
    x.append(i[1][1])
    m.append(i[1][2])
y = [i for i in range(len(x))]
plt.plot(y, x, "--", c='r')
plt.plot(y, m, "--", c='b')

mesh_min = mesh[1:-1, 1:-1].min() - 0.01
mesh_max = mesh[1:-1, 1:-1].max()
plt.imshow(mesh, vmin=mesh_min, vmax=mesh_max)
