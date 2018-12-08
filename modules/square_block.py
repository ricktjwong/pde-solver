#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 00:47:27 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

h = 1E-3                    # Step size h (in m)
k = 150                     # Conductivity of silicon Microchip in W/m K
T_a = 20                    # Ambient temperature
alpha = h * 2.62 / k        # Constant for natural convection
del_T = 0.25 * h ** 2 * 500 * 1E6 / k
s = 1                      # Introduce a scaling factor
rows = 10 * s
cols = 10 * s


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
mesh[1:-1, 1:-1] = 20      # Initialise the inner mesh with T > T_a
mesh = update_boundaries(mesh).copy()

all_mesh = []

n = 0
while (n < 100000):
    update = mesh.copy()
    for j in range(1, cols-1):
        for i in range(1, rows-1):
            update[i][j] = 1/4 * (mesh[i-1][j] + mesh[i+1][j] 
                                + mesh[i][j-1] + mesh[i][j+1])
    update[4*s:6*s, 4*s:6*s] += del_T
#    if np.linalg.norm(update[1:-1, 1:-1]/mesh[1:-1, 1:-1] - 1) < 1E-4:
#        break
    update = update_boundaries(update).copy()
    mesh = update.copy()
    all_mesh.append(mesh)
    n += 1

x = []
m = []
for i in all_mesh:
    x.append(i[1][1])
    m.append(i[1][2])
y = [i for i in range(len(x))]
plt.figure(2)
plt.plot(y, x, "--", c='r')
plt.plot(y, m, "--", c='b')

T_min = mesh[1:-1, 1:-1].min() - 0.01
T_max = mesh[1:-1, 1:-1].max()
masked = np.ma.masked_where(mesh < 0.01, mesh)
ax = plt.subplot(111)
im = ax.imshow(masked, cmap='hot', vmin=T_min, vmax=T_max)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.xlabel('$T$ / K', labelpad=20)
plt.colorbar(im, cax=cax)

#masked = np.ma.masked_where(mesh < 0.01, mesh)
#plt.imshow(masked, cmap='rainbow')
