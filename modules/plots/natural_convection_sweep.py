#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 18:13:38 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 30
plt.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['figure.figsize'] = (10, 7)

vary_nfins = np.loadtxt("../../data/natural_convection/RB_SOR/vary_nfins.txt")
x = vary_nfins[:,0]
T = vary_nfins[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=200)
#plt.savefig("vary_nfins.pdf", dpi=3000)

vary_height = np.loadtxt("../../data/natural_convection/RB_SOR/vary_fin_height.txt")
x = vary_height[:,0]
T = vary_height[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=200)
#plt.savefig("vary_fin_height.pdf", dpi=3000)

vary_thickness = np.loadtxt("../../data/natural_convection/RB_SOR/vary_fin_thickness.txt")
x = vary_thickness[:,0]
T = vary_thickness[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=200)
#plt.savefig("vary_fin_thickness.pdf", dpi=3000)

vary_sep = np.loadtxt("../../data/natural_convection/RB_SOR/vary_fin_sep.txt")
x = vary_sep[:,0]
T = vary_sep[:,1]
plt.figure()
plt.scatter(x, T, marker="x", c='r', s=200)
#plt.savefig("vary_fin_sep.pdf", dpi=3000)
