#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:13:10 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import modules.utils.solvers as solv
import modules.ceramic_cooling as cc
import modules.heat_structure as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 24
plt.rcParams['lines.linewidth'] = 1.4
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

"""
Question 3:
No Heat Sink
Ceramic cooling

Estimated run time is 70s for the Red-Black SOR method
Converging from ambient temperature, we get 7095.675215837743, with 218,998
iterations and when converging from the top (i.e. setting initial temperature
of the structure to 9000 degrees, we get 7100.049272489449). The error is thus
calculated to be half the difference of the two.
"""

#mp = cc.Microprocessor(5, conv_ratio = 1E-8, solver = solv.gauss_seidel)
#start = time.time()
#final_temp, n = mp.solve_mesh()
#end = time.time()
#print(final_temp, n)
#print(end - start)

#T = np.load("all_temps.npy")
#x = [i for i in range(len(T))]

#plt.figure()
#plt.plot(x, T, "--", c='r')

#plt.figure()
#final_mesh = np.zeros((mp.rows, mp.cols))
#mp.update_nonboundaries(all_mesh[-1], final_mesh)
#
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#ax = plt.subplot(111)
#im = ax.imshow(masked, cmap='rainbow')
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.xlabel('$T$ / K', labelpad=20)
#plt.colorbar(im, cax=cax)


"""
Natural convection
b: fin separation
c: fin thickness
f_h: fin height
scale, b, c, f_h, n_fins
"""

#hs = hst.HeatStructure(2, 5, 2, 30, 5, convection_type="natural")
#start = time.time()
#all_mesh = hs.solve_mesh()
#end = time.time()
#print(end - start)
#print(np.mean(all_mesh[-1][hs.m_idx_y1:hs.m_idx_y2, hs.m_idx_x1:hs.m_idx_x2]))

"""
Force convection
"""

#hs = hst.HeatStructure(2, 5, 1, 30, 5, conv_ratio=1E-6,
hs = hst.HeatStructure(4, 1, 1, 100, 100, conv_ratio=1E-8,
                       convection_type="forced", solver=solv.red_black_SOR)
start = time.time()
temps, n = hs.solve_mesh()
end = time.time()
print(end - start)
print(temps[-1], n)

plt.figure()
final_mesh = np.zeros((hs.rows, hs.cols))
hs.update_nonboundaries(hs.mesh, final_mesh)

masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
ax = plt.subplot(111)
im = ax.imshow(masked, cmap='rainbow')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.xlabel('$T$ / $^o$C', labelpad=20)
plt.colorbar(im, cax=cax)

#plt.figure(2)
#x = [i for i in range(1, len(temps)+1)]
#plt.plot(x, temps, "--", c='r')
#
#plt.figure(3)
#final_mesh = np.zeros((hs.rows, hs.cols))
#hs.update_nonboundaries(hs.mesh, final_mesh)
#
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#plt.imshow(masked, cmap="rainbow")
#plt.show()
