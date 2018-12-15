#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:13:10 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import utils.solvers as solv
import modules.ceramic_cooling as cc
import modules.heat_structure as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
plt.rcParams['lines.linewidth'] = 1.4
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

"""
Question 3:
No Heat Sink
Ceramic cooling
Estimated run time is 10 minutes
"""

#mp = cc.Microprocessor(5, conv_ratio = 1E-8, solver = solv.jacobi_solver)
#all_mesh, temps, final_temp, n = mp.solve_mesh()
#print(final_temp)
#np.save("average_temps", temps)
#np.save("all_mesh", all_mesh)
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

#final = []
#fins = np.arange(15, 50, 5)
#for i in range(12, 13):
#hs = hst.HeatStructure(1, 2, 2, 30, 40, conv_ratio=1E-6,
#                       convection_type="forced", solver=solv.jacobi_solver)
#start = time.time()
#final_temp, n = hs.solve_mesh()
#all_mesh = hs.solve_mesh()    
#end = time.time()
#print(end - start)
#    final.append([i, final_temp, n, (end - start)])
#print(final_temp, n)
#np.savetxt("run1.txt", final)
#print(np.mean(all_mesh[-1][hs.m_idx_y1:hs.m_idx_y2, hs.m_idx_x1:hs.m_idx_x2]))

#x = []
#for i in all_mesh:
#    x.append(np.mean(i[hs.m_idx_y1:hs.m_idx_y2, hs.m_idx_x1:hs.m_idx_x2]))
#y = [i for i in range(len(x))]
#
#print(x[-1])
#
#plt.figure(2)
#plt.plot(y, x, "--", c='r')
#
#plt.figure(3)
#final_mesh = np.zeros((hs.rows, hs.cols))
#hs.update_nonboundaries(hs.mesh, final_mesh)
#
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#plt.imshow(masked, cmap="rainbow")
#plt.show()
