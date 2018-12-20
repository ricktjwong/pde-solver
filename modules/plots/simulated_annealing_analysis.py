#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 15:02:35 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt


costs = np.load("../../data/sim_annealing/costs_forced.npy")
actions = np.load("../../data/sim_annealing/actions_forced.npy")

b, c, f, n = actions[:,0], actions[:,1], actions[:,2], actions[:,3]
x = np.arange(len(b))

fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True)
ax = axs[0]
ax.set_xlabel('Number of accepted values')
ax.set_ylabel('Parameter unit')    
ax.scatter(x, b, marker="x", c='r', label='Fin separation b')
ax.scatter(x, c, marker="x", c='b', label='Fin thickness c')
ax.scatter(x, f, marker="x", c='g', label='Fin height a')
ax.scatter(x, n, marker="x", c='y', label='Number of fins n')
ax.legend()
ax.set_title('Parameter Change for Every Accepted Iteration')

ax = axs[1]
ax.set_xlabel('Number of accepted values')
ax.set_ylabel('Cost value ($^o$C)')
ax.scatter(x, costs, marker="x", c='black')
ax.set_title('Cost for Every Accepted Iteration')

fig.suptitle('Simulated Annealing')
