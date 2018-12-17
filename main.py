#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:13:10 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import modules.automated_runs.natural_convection_runs as ncr
import modules.automated_runs.forced_convection_runs as fcr
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
Ceramic-Microprocessor cooling

Estimated run time is 600s for the Jacobi iterative method
Estimated run time is 70s for the Red-Black SOR method
"""

print("###############################################")
print("Question 3: Average operating temperature of microchip under" +
      "natural convection")
print("###############################################")
      
mp = cc.Microprocessor(5, conv_ratio=1E-8, solver=solv.jacobi_solver)
T, n = mp.solve_mesh()
print("The average microprocessor temperature under natural convection and" + \
      " without a heat sink is: {0}.".format(T[-1]))

plt.figure()
x = [i for i in range(len(T))]
plt.plot(x, T, "--", c='r')

plt.figure()
final_mesh = np.zeros((mp.rows, mp.cols))
mp.update_nonboundaries(mp.mesh, final_mesh)
masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
ax = plt.subplot(111)
im = ax.imshow(masked, cmap='hot')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.xlabel('$T$ / $^o$C', labelpad=20)
plt.colorbar(im, cax=cax)

"""
Question 4:
Heat Sink with Natural Convection
Estimated total run time is 6 hours
"""

print("###############################################")
print("Question 4: Vary params to see which works the best for decreasing" +
      "the average temperature of the microprocessor under natural convection")
print("###############################################")

# Vary b, c, f_h, n_fins
# Sweeping all four params takes about 3 hours
ncr.sweep_params()
# Vary only height and number of fins. THis will take about 3 hours
ncr.minimise_temp()


"""
Question 5
Heat Sink with Forced Convection
Estimated total run time is 2 hours
"""

print("###############################################")
print("Question 5: Comparison of the microprocessor average temperature under"
      + "natural and forced convection using the Red-Black SOR method")
print("###############################################")

# This will take about 2 minutes to run
hs = hst.HeatStructure(2, 5, 1, 30, 5, conv_ratio=1E-6,
                       convection_type="natural", solver=solv.red_black_SOR)
T, n = hs.solve_mesh()
print("The average microprocessor temperature under natural convection is: " +
      "{0}.".format(T[-1]))

plt.figure()
final_mesh = np.zeros((hs.rows, hs.cols))
hs.update_nonboundaries(hs.mesh, final_mesh)

masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
ax = plt.subplot(111)
im = ax.imshow(masked, cmap='hot')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.xlabel('$T$ / $^o$C', labelpad=20)
plt.colorbar(im, cax=cax)

hs = hst.HeatStructure(2, 5, 1, 30, 5, conv_ratio=1E-6,
                       convection_type="forced", solver=solv.red_black_SOR)
T, n = hs.solve_mesh()
print("The average microprocessor temperature under forced convection is: " +
      "{0}.".format(T[-1]))

plt.figure()
final_mesh = np.zeros((hs.rows, hs.cols))
hs.update_nonboundaries(hs.mesh, final_mesh)

masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
ax = plt.subplot(111)
im = ax.imshow(masked, cmap='hot')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.xlabel('$T$ / $^o$C', labelpad=20)
plt.colorbar(im, cax=cax)

# This will take about 1 hour to run
# Do a sweep to increase height and number of fins to get T below 80 degrees
fcr.sweep_params()
