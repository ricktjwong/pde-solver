#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 01:01:31 2018

@author: ricktjwong
"""

import numpy as np
import time
import modules.heat_structure as hst

"""
Natural convection
scale, b, c, f_h, n_fins
"""

# Vary n_fins
#final = []
#n_fins = np.arange(5, 65, 5)
#for i in n_fins:
#    hs = hst.HeatStructure(2, 5, 1, 30, i, convection_type="natural")
#    start = time.time()
#    final_temp, n = hs.solve_mesh()
#    end = time.time()
#    final.append([i, final_temp, n, (end - start)])
#np.savetxt("natural_convection_vary_nfins.txt", final)
#
## Vary fin thickness c
#final = []
#fin_thickness = np.arange(1, 21, 1)
#for i in fin_thickness:
#    hs = hst.HeatStructure(2, 5, i, 30, 5, convection_type="natural")
#    start = time.time()
#    final_temp, n = hs.solve_mesh()
#    end = time.time()
#    final.append([i, final_temp, n, (end - start)])
#np.savetxt("natural_convection_vary_fin_thickness.txt", final)
#
## Vary fin height f_h
#final = []
#fin_heights = np.arange(30, 145, 5)
#for i in fin_heights:
#    hs = hst.HeatStructure(2, 5, 1, i, 5, convection_type="natural")
#    start = time.time()
#    final_temp, n = hs.solve_mesh()
#    end = time.time()
#    final.append([i, final_temp, n, (end - start)])
#    print(final)
#np.savetxt("natural_convection_vary_fin_height.txt", final)

# Vary fin separation b
final = []
fin_seps = np.arange(1, 26, 1)
for i in fin_seps:
    hs = hst.HeatStructure(2, i, 1, 30, 12, convection_type="natural")
    start = time.time()
    final_temp, n = hs.solve_mesh()
    end = time.time()
    final.append([i, final_temp, n, (end - start)])
    print(final)    
np.savetxt("natural_convection_vary_fin_sep.txt", final)
