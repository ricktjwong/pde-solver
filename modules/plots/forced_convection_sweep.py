#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:13:38 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_sweeps():
    """
    Forced convection
    scale, b, c, f_h, n_fins
    """
    # scale, b, c, f_h = 4, 1, 1, 100
    vary_nfins = np.loadtxt("data/forced_convection/vary_nfins.txt")
    x = vary_nfins[:,0]
    T = vary_nfins[:,1]
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True)
    ax = axs[0]
    ax.set_xlabel('Number of fins')
    ax.set_ylabel('Average Microprocessor Temperature ($^o$C)')
    ax.set_title('Average Microprocessor Temperature with Varying Number of Fins')
    ax.scatter(x, T, marker="x", c='r', s=200,
                label="Temperature at each fin number")
    ax.legend()
    
    # scale, b, c, n_fins = 4, 1, 1, 100
    vary_height = np.loadtxt("data/forced_convection/vary_f_h.txt")
    x = vary_height[:,0]
    T = vary_height[:,1]
    ax = axs[1]
    ax.set_xlabel('Fin height (mm)')
    ax.set_ylabel('Average Microprocessor Temperature ($^o$C)')
    ax.set_title('Average Microprocessor Temperature with Varying Fin Height')
    ax.scatter(x, T, marker="x", c='r', s=200,
                label="Temperature at each fin height")
    ax.legend()

    fig.suptitle('Forced convection sweep')
