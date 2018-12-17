#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:13:38 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_sweeps():
    vary_nfins = np.loadtxt("data/natural_convection/RB_SOR/vary_nfins.txt")
    x = vary_nfins[:,0]
    T = vary_nfins[:,1]
    fig, axs = plt.subplots(nrows=2, ncols=3)
    ax = axs[0, 0]
    ax.set_xlabel("Number of fins")  
    ax.set_ylabel("Average Microprocessor Temperature ($^o$C)")
    ax.scatter(x, T, marker="x", c='r', s=200)
    
    vary_height = np.loadtxt("data/natural_convection/RB_SOR/vary_fin_height.txt")
    x = vary_height[:,0]
    T = vary_height[:,1]
    ax = axs[0, 1]
    ax.set_xlabel("Fin height (mm)")
    ax.set_ylabel("Average Microprocessor Temperature ($^o$C)")
    ax.scatter(x, T, marker="x", c='r', s=200)
    
    vary_thickness = np.loadtxt("data/natural_convection/RB_SOR/vary_fin_thickness.txt")
    x = vary_thickness[:,0]
    T = vary_thickness[:,1]
    ax = axs[0, 2]
    ax.set_xlabel("Fin thickness (mm)")
    ax.set_ylabel("Average Microprocessor Temperature ($^o$C)")    
    ax.scatter(x, T, marker="x", c='r', s=200)
    
    vary_sep = np.loadtxt("data/natural_convection/RB_SOR/vary_fin_sep.txt")
    x = vary_sep[:,0]
    T = vary_sep[:,1]
    ax = axs[1, 0]
    ax.set_xlabel("Fin Separation (mm)")
    ax.set_ylabel("Average Microprocessor Temperature ($^o$C)")    
    ax.scatter(x, T, marker="x", c='r', s=200)
    
    # Sweep f_h from 30 to 200
    # Scale, b, c, n_fins = 2, 1, 1, 50
    vary_height = np.loadtxt("data/natural_convection/vary_f_h.txt")
    x = vary_height[:,0]
    T = vary_height[:,1]
    ax = axs[1, 1]
    ax.set_xlabel("Fin height (mm)")
    ax.set_ylabel("Average Microprocessor Temperature ($^o$C)")    
    ax.scatter(x, T, marker="x", c='r', s=200)
    
    # Sweep n_fins from 30 to 200
    # Scale, b, c, n_fins = 2, 1, 1, 50
    vary_nfins = np.loadtxt("data/natural_convection/vary_nfins.txt")
    x = vary_nfins[:,0]
    T = vary_nfins[:,1]
    ax = axs[1, 2]
    ax.set_xlabel("Number of fins")
    ax.set_ylabel("Average Microprocessor Temperature ($^o$C)")    
    ax.scatter(x, T, marker="x", c='r', s=200)

    fig.suptitle('Natural convection sweep')
