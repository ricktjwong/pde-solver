#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:45:07 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
plt.rcParams['lines.linewidth'] = 1.4
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

T = np.loadtxt("data/ceramic_cooling/average_temp.txt")
x = np.arange(1, len(T)+1, 1)
#plt.figure()
#plt.plot(x, T, c='r')
#plt.savefig("ceramic_cooling_T.pdf", dpi=3000)

#final_mesh = np.loadtxt("data/ceramic_cooling/final_mesh.txt")
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#plt.figure()
#ax = plt.subplot(111)
#im = ax.imshow(masked, cmap='rainbow')
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.xlabel('$T$ / K', labelpad=20)
#plt.colorbar(im, cax=cax)