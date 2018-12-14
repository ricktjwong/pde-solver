#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 03:55:10 2018

@author: ricktjwong
"""

import numpy as np
import time
import modules.heat_structure as hst

"""
Forced convection
scale, b, c, f_h, n_fins
"""

# Vary the scale and test convergence temperature
final = []
scale = np.arange(2, 16, 1)
for i in scale:
    if 1E-6/(10 ** i) > 1E-14:
        ratio = 1E-6/(10 ** i)
    else:
        ratio = 1E-14
    hs = hst.HeatStructure(i, 1, 1, 30, 30, convection_type="forced",
                           conv_ratio=ratio)
    start = time.time()
    final_temp, n = hs.solve_mesh()
    end = time.time()
    final.append([i, final_temp, n, (end - start)])
    print(final)
np.savetxt("forced_convection_vary_scale.txt", final)
