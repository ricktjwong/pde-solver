#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:13:10 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import modules.heat_structure as hst

"""
Natural convection
scale, b, c, f_h, n_fins
"""

#hs = hst.HeatStructure(2, 5, 2, 30, 5, convection_type="natural")
#start = time.time()
#all_mesh = hs.solve_mesh()
#end = time.time()
#print(end - start)
#print(np.mean(all_mesh[-1][hs.m_idx_y1:hs.m_idx_y2, hs.m_idx_x1:hs.m_idx_x2]))

"""
Force convection
"""

final = []
for i in range(1, 2):
    hs = hst.HeatStructure(i, 5, 2, 30, 5, convection_type="force")
    start = time.time()
    final_temp, n = hs.solve_mesh()
    end = time.time()
    print(end - start)
    final.append([i, final_temp, n])
    print(final)
np.savetxt("run1.txt", final)
#print(np.mean(all_mesh[-1][hs.m_idx_y1:hs.m_idx_y2, hs.m_idx_x1:hs.m_idx_x2]))

#x = []
#for i in all_mesh:
#    x.append(np.mean(i[hs.m_idx_y1:hs.m_idx_y2, hs.m_idx_x1:hs.m_idx_x2]))
#y = [i for i in range(len(x))]

#plt.figure(2)
#plt.plot(y, x, "--", c='r')
#
#plt.figure(3)
#final_mesh = np.zeros((hs.rows, hs.cols))
#hs.update_nonboundaries(all_mesh[-1], final_mesh)
#
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#plt.imshow(masked, cmap="rainbow")
#plt.show()
