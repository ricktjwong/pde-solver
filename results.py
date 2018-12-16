#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 01:40:21 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import modules.utils.solvers as solv
import modules.ceramic_cooling as cc
import modules.heat_structure as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.rcParams['lines.linewidth'] = 1.4
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

"""
Question 3:
No Heat Sink
Ceramic cooling
"""

#mp = cc.Microprocessor(5, conv_ratio = 1E-8, solver = solv.jacobi_solver)
#T = np.load("data/ceramic_cooling/average_temps.npy")
#x = np.arange(1, len(T)+1, 1)
#plt.figure(figsize=(10, 7))
#plt.plot(x, T, c='r', linestyle='--')
#plt.savefig("ceramic_cooling_T.pdf", dpi=3000)

#all_mesh = np.load("data/ceramic_cooling/all_mesh.npy")

#plt.figure()
#final_mesh = np.zeros((mp.rows, mp.cols))
#mp.update_nonboundaries(all_mesh[100], final_mesh)
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#ax = plt.subplot(111)
#im = ax.imshow(masked, cmap='hot')
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.xlabel('$T$ / K', labelpad=20)
#plt.colorbar(im, cax=cax)
