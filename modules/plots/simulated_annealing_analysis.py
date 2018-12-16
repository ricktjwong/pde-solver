#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 15:02:35 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

costs = np.load("../../data/costs0.00025.npy")
actions = np.load("../../data/action0.00025.npy")

b, c, f, n = actions[:,0], actions[:,1], actions[:,2], actions[:,3]
x = np.arange(len(b))
plt.scatter(x, b, marker="x", c='r')
plt.scatter(x, c, marker="x", c='b')
plt.scatter(x, f, marker="x", c='g')
plt.scatter(x, n, marker="x", c='y')

plt.figure()
plt.scatter(x, costs, marker="x", c='black')
