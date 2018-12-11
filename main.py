#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:13:10 2018

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import modules.natural_convection as nc
import modules.force_convection as fc

#mesh = fc.setup()
#c_mesh, m_mesh, fb_mesh, f_mesh = fc.update_all_boundaries(mesh)
#mesh = fc.update_mesh(mesh, c_mesh, m_mesh, fb_mesh, f_mesh).copy()
#
#start = time.time()
#all_mesh = fc.solve_mesh(mesh, 1E-6)
#end = time.time()
#print(end - start)
#print(np.mean(all_mesh[-1][fc.m_idx_y1:fc.m_idx_y2, fc.m_idx_x1:fc.m_idx_x2]))

mesh = nc.setup()
c_mesh, m_mesh, fb_mesh, f_mesh = nc.update_all_boundaries(mesh)
mesh = nc.update_mesh(mesh, c_mesh, m_mesh, fb_mesh, f_mesh).copy()

start = time.time()
all_mesh = nc.solve_mesh(mesh, 1E-6)
end = time.time()
print(end - start)
print(np.mean(all_mesh[-1][nc.m_idx_y1:nc.m_idx_y2, nc.m_idx_x1:nc.m_idx_x2]))

#x = []
#for i in all_mesh:
#    x.append(np.mean(i[fc.m_idx_y1:fc.m_idx_y2, fc.m_idx_x1:fc.m_idx_x2]))
#y = [i for i in range(len(x))]
#
#plt.figure(2)
#plt.plot(y, x, "--", c='r')

#plt.figure(3)
#final_mesh = np.zeros((fc.rows, fc.cols))
#fc.update_nonboundaries(all_mesh[-1], final_mesh)

#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#plt.imshow(masked, cmap="rainbow")
#plt.show()
