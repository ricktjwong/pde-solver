#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 01:01:31 2018

@author: ricktjwong
"""

import sys
sys.path.append("../../")
import numpy as np
import time
import modules.utils.solvers as solv
import modules.heat_structure as hst

"""
Natural convection
scale, b, c, f_h, n_fins
"""

def sweep_params():
    """
    Sweep all four parameters, keeping the other three constant
    """
    print("hi")
    # Vary n_fins
    final = []
    n_fins = np.arange(5, 65, 5)
    for i in n_fins:
        hs = hst.HeatStructure(2, 5, 1, 30, i, convection_type="natural",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
    np.savetxt("natural_convection_vary_nfins.txt", final)
    
    # Vary fin thickness c
    final = []
    fin_thickness = np.arange(1, 31, 1)
    for i in fin_thickness:
        hs = hst.HeatStructure(2, 5, i, 30, 5, convection_type="natural",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
    np.savetxt("natural_convection_vary_fin_thickness.txt", final)
    
    # Vary fin height f_h
    final = []
    fin_heights = np.arange(30, 205, 5)
    for i in fin_heights:
        hs = hst.HeatStructure(2, 5, 1, i, 5, convection_type="natural",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
    np.savetxt("natural_convection_vary_fin_height.txt", final)
    
    # Vary fin separation b
    final = []
    fin_seps = np.arange(5, 41, 1)
    for i in fin_seps:
        hs = hst.HeatStructure(2, i, 1, 30, 5, convection_type="natural",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
    np.savetxt("natural_convection_vary_fin_sep.txt", final)

"""
Sweep to try to get a temperature of below 100 degrees
"""

def minimise_temp():
    """
    Sweep height and fins to try to get below 100 degrees
    """
    # Sweep f_h
    final = []
    f_h = np.arange(200, 205, 5)
    for i in f_h:
        hs = hst.HeatStructure(2, 1, 1, i, 50, convection_type="natural",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
        print(final)
    np.savetxt("natural_runs_3.txt", final)
        
    # Sweep n_fins
    final = []
    nfins = np.arange(30, 205, 5)
    for i in nfins:
        hs = hst.HeatStructure(2, 1, 1, 50, i, convection_type="natural",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
        print(final)
    np.savetxt("natural_runs_vary_nfins.txt", final)
