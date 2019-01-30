#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 13:16:18 2019

@author: ricktjwong
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import modules.utils.solvers as solv
import modules.heat_structure as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['figure.figsize'] = (10, 10)

fig = plt.figure()
ax = plt.subplot(111)
a = np.random.random((76, 52))
final_mesh = np.zeros((76, 52))
masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
im = plt.imshow(a, interpolation='nearest', cmap='hot', vmin=20,
                vmax=100, norm=matplotlib.colors.LogNorm(clip=True))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.xlabel('$T$ / $^o$C', labelpad=20)
plt.colorbar(im, cax=cax)

hs = hst.HeatStructure(2, 5, 1, 30, 5, im, conv_ratio=1E-6,
                       convection_type="forced", solver=solv.SOR_solver)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, hs.animate, init_func=hs.initialise_animation,
                               frames=1000, interval=100, blit=False)

anim.save('basic_animation_hot.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
