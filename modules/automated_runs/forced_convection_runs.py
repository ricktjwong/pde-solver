  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 19:08:36 2018

@author: ricktjwong
"""

import sys
sys.path.append("../../")
import numpy as np
import time
import modules.utils.solvers as solv
import modules.heat_structure as hst


def sweep_params():
    """
    Find minimum temperature by performing a sweep of height, fins and wind speed
    """
    # Vary fin height
    final = []
    f_h = np.arange(100, 150, 5)
    for i in f_h:
        hs = hst.HeatStructure(4, 1, 1, i, 100, convection_type="forced",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
    np.savetxt("forced_convection_vary_f_h.txt", final)
    
    # Vary n_fins
    final = []
    n_fins = np.arange(100, 150, 5)
    for i in n_fins:
        hs = hst.HeatStructure(4, 1, 1, 100, i, convection_type="forced",
                               solver=solv.red_black_SOR)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
    np.savetxt("forced_convection_vary_nfins.txt", final)
    
    # Vary wind_speed
    final = []
    v = np.arange(20, 31, 1)
    for i in v:
        hs = hst.HeatStructure(4, 1, 1, 100, 100, convection_type="forced",
                               solver=solv.red_black_SOR, wind_speed = i)
        start = time.time()
        final_temp, n = hs.solve_mesh()
        end = time.time()
        final.append([i, final_temp, n, (end - start)])
        print(final)
    np.savetxt("forced_convection_vary_windspeed.txt", final)
