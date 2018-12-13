#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:13:38 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

vary_nfins = np.loadtxt("./data/natural_convection/vary_nfins.txt")
x = vary_nfins[:,0]
T = vary_nfins[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r')

vary_height = np.loadtxt("./data/natural_convection/vary_fin_height.txt")
x = vary_height[:,0]
T = vary_height[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r')

vary_thickness = np.loadtxt("./data/natural_convection/vary_fin_thickness.txt")
x = vary_thickness[:,0]
T = vary_thickness[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r')

vary_sep = np.loadtxt("./data/natural_convection/vary_fin_sep.txt")
x = vary_sep[:,0]
T = vary_sep[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r')
