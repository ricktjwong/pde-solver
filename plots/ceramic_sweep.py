#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:13:38 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['figure.figsize'] = (10, 7)

vary_nfins = np.loadtxt("../data/natural_convection/vary_nfins.txt")
x = vary_nfins[:,0]
T = vary_nfins[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=100)
#plt.savefig("vary_nfins.pdf", dpi=3000)

vary_height = np.loadtxt("../data/natural_convection/vary_fin_height.txt")
x = vary_height[:,0]
T = vary_height[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=100)
#plt.savefig("vary_fin_height.pdf", dpi=3000)

vary_thickness = np.loadtxt("../data/natural_convection/vary_fin_thickness.txt")
x = vary_thickness[:,0]
T = vary_thickness[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=100)

vary_sep = np.loadtxt("../data/natural_convection/vary_fin_sep.txt")
x = vary_sep[:,0]
T = vary_sep[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=100)

x = np.array([[2, 82.66034679658038, 90301, 129.9484899044037],
             [3, 93.12628579331323, 301297, 731.2707726955414],
             [4, 98.65914099965012, 697794, 1913.542496919632],
             [5, 101.97653010470077, 1296485, 10307.263709783554]])
s = x[:,0]
T = x[:,1]
plt.figure()
plt.scatter(s, T, marker="x", c='r', s=100)
