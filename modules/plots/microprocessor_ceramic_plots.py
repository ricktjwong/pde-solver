#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:45:07 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import modules.utils.solvers as solv
import modules.ceramic_cooling as cc
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_figures():
    mp = cc.Microprocessor(5, conv_ratio = 1E-8, solver = solv.jacobi_solver)
    T = np.load("data/ceramic_microprocessor/average_temps.npy")
    x = np.arange(1, len(T)+1, 1)
    plt.figure(figsize=(10, 7))
    plt.xlabel('Number of iterations')
    plt.ylabel('Average Microprocessor Temperature ($^o$C)')
    plt.title('Average Microprocessor Temperature with Iteration')
    plt.plot(x, T, c='r', linestyle='--')
    
    all_mesh = np.load("data/ceramic_microprocessor/all_mesh.npy")
    
    plt.figure()
    final_mesh = np.zeros((mp.rows, mp.cols))
    mp.update_nonboundaries(all_mesh[-1], final_mesh)
    masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
    ax = plt.subplot(111)
    ax.set_xlabel("x-position (grid counts")
    ax.set_ylabel("y-position (grid counts")
    im = ax.imshow(masked, cmap='hot')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.xlabel('$T$ / $^o$C', labelpad=20)
    plt.colorbar(im, cax=cax)
