#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 03:09:37 2018

@author: ricktjwong
"""

import numpy as np
from utils.jacobi_method import jacobi_solver
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

k_m = 150                   # Conductivity of silicon Microchip in W/m K
k_c = 230                   # Conductivity of ceramic block in W/m K
k_a = 248                   # Conductivity of aluminium fin in W/m K
T_a = 20                    # Ambient temperature

class Microprocessor():
    def __init__(self, scale, conv_ratio = 1E-6):
        self.scale, self.conv = scale, conv_ratio
        # Step size h (in m)
        self.h = 1E-3 / self.scale
        # Change in temperature of microprocessor every s
        self.rho = 0.25 * self.h ** 2 * 500 * 1E6 / k_m
        self.phi_m = self.h * 2.62 / k_m 
        self.phi_c = self.h * 2.62 / k_c
        self.rows = 3 * self.scale + 2
        self.cols = 20 * self.scale + 2
        self.mesh = np.zeros((self.rows, self.cols)) 
        self.initialise_indices()
        self.initialise_temp()
        self.initialise_boundaries()
    
    def initialise_indices(self):
        """
        Initialise indices
        """        
        # Initialise row and column indices for ceramic block
        self.c_idx_x1, self.c_idx_x2 = 1, self.cols - 1
        self.c_idx_y1, self.c_idx_y2 = 1, 2 * self.scale + 1
        # Initialise row and column indices for microprocessor block
        self.m_idx_x1 = int((self.cols - 14 * self.scale)/2)
        self.m_idx_x2 = self.m_idx_x1 + 14 * self.scale
        self.m_idx_y1, self.m_idx_y2 = \
        self.c_idx_y2, self.c_idx_y2 + self.scale

    def initialise_temp(self):
        """
        Set start temperature of the components to ambient temperature T_a
        """
        # Initialise the ceramic block at T_a
        self.mesh[self.c_idx_y1:self.c_idx_y2,
                  self.c_idx_x1:self.c_idx_x2] = T_a
        # Initialise the microprocessor at T_a
        self.mesh[self.m_idx_y1:self.m_idx_y2,
                  self.m_idx_x1:self.m_idx_x2] = T_a
        
    def initialise_boundaries(self):
        c_mesh, m_mesh = self.update_all_boundaries(self.mesh)
        self.update_mesh(self.mesh, c_mesh, m_mesh)
    
    def update_boundaries(self, m, c):
        """
        Calculate the values for boundaries using the Central Difference Scheme
        (CDS). Each mesh (m) represents a rectangular component of the entire
        structure
        """
        rows = m.shape[0]
        cols = m.shape[1]
        # Extract the non-boundary mesh points for four sides of the block
        l_mesh = m[1:-1, 0:3]
        r_mesh = m[1:-1, cols-3:]
        u_mesh = m[0:3, 1:-1]
        d_mesh = m[rows-3:, 1:-1]
        # Update the boundary points using the CDS for four sides of block
        l_mesh[:,0] = l_mesh[:,-1] - c * (l_mesh[:,1] - T_a) ** (4/3)
        r_mesh[:,-1] = r_mesh[:,0] - c * (r_mesh[:,1] - T_a) ** (4/3)
        u_mesh[0] = u_mesh[-1] - c * (u_mesh[1] - T_a) ** (4/3)
        d_mesh[-1] = d_mesh[0] - c * (d_mesh[1] - T_a) ** (4/3)
        # Update the original mesh with the new boundary values
        m[1:-1, 0:3] = l_mesh
        m[1:-1, cols-3:] = r_mesh
        m[0:3, 1:-1] = u_mesh
        m[rows-3:, 1:-1] = d_mesh
    
    def update_all_boundaries(self, m):
        """
        Update boundaries for each building block of the entire structure
        """
        c_mesh = m[self.c_idx_y1-1:self.c_idx_y2+1,
                   self.c_idx_x1-1:self.c_idx_x2+1].copy()
        self.update_boundaries(c_mesh, self.phi_c)
        m_mesh = m[self.m_idx_y1-2:self.m_idx_y2+1,
                   self.m_idx_x1-1:self.m_idx_x2+1].copy()
        self.update_boundaries(m_mesh, self.phi_m)
        return c_mesh, m_mesh
    
    def update_mesh(self, m, c_mesh, m_mesh):
        """
        Update original mesh with the boundary values of different components. 
        The components are stacked in an order that preserves the correct
        boundaries. In this case we stack in the order: ceramic, microchip
        """
        m[self.c_idx_y1 - 1:self.c_idx_y2 + 1,
          self.c_idx_x1 - 1:self.c_idx_x2 + 1] = c_mesh.copy()
        m[self.m_idx_y1:self.m_idx_y2 + 1,
          self.m_idx_x1-1:self.m_idx_x2 + 1] = m_mesh[2:].copy()
        # Handle the overlap boundaries
        m[self.m_idx_y1][self.m_idx_x1 - 1] = \
        (m[self.m_idx_y1][self.m_idx_x1 - 1] + \
         c_mesh[-1][self.m_idx_x1 - self.c_idx_x1]) / 2
        m[self.m_idx_y1][self.m_idx_x2] = \
        (m[self.m_idx_y1][self.m_idx_x2] + \
         c_mesh[-1][self.m_idx_x2 - self.c_idx_x2 - 1]) / 2        
        return m
    
    def update_nonboundaries(self, m, n):
        """
        Update all the mesh points for the structure excluding the boundary
        values. m is the original mesh, n is the new mesh to be returned
        """
        # Update matrix for the ceramic block    
        n[self.c_idx_y1:self.c_idx_y2, self.c_idx_x1:self.c_idx_x2] = \
        m[self.c_idx_y1:self.c_idx_y2, self.c_idx_x1:self.c_idx_x2]
        # Update matrix for the microprocessor
        n[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2] = \
        m[self.m_idx_y1:self.m_idx_y2, self.m_idx_x1:self.m_idx_x2]
        
    def jacobi_solver(self, x):
        x_u = np.roll(x, -1, 0)
        x_u[-1] = np.zeros((1, x.shape[1]))
        x_d = np.roll(x, 1, 0)
        x_d[0] = np.zeros((1, x.shape[1]))        
        x_r = np.roll(x, 1, 1)
        x_r[:,0] = np.zeros((1, x.shape[0]))        
        x_l = np.roll(x, -1, 1)
        x_l[:,-1] = np.zeros((1, x.shape[0]))        
        y = (x_d + x_u + x_r + x_l) / 4
        return y    
    
    def solve_mesh(self):
        """
        Iterate and solve for the temperature at every mesh point till the 
        change in average temperature of the microprocessor is below the
        convergence_ratio
        """
        n = 0
        all_mesh = []
        while (True):
            update = self.mesh.copy()
            out = jacobi_solver(update)
            self.update_nonboundaries(out, update)
            update[self.m_idx_y1:self.m_idx_y2,
                   self.m_idx_x1:self.m_idx_x2] += self.rho
            m_mean_temp_1 = np.mean(self.mesh[self.m_idx_y1:self.m_idx_y2,
                                              self.m_idx_x1:self.m_idx_x2])
            norm_temp_1 = np.linalg.norm(self.mesh)
            norm_temp_2 = np.linalg.norm(update)
            if norm_temp_2 / norm_temp_1 - 1 < self.conv: break
            c_mesh, m_mesh = self.update_all_boundaries(update)
            update = self.update_mesh(update, c_mesh, m_mesh).copy()
            self.mesh = update.copy()
            if n % 10000 == 0:
                print(n)
            all_mesh.append(self.mesh)
            n += 1
        return all_mesh
#        return m_mean_temp_1, n

mp = Microprocessor(5, conv_ratio = 1E-8)
#plt.imshow(mp.mesh)
all_mesh = mp.solve_mesh()
print(np.mean(all_mesh[-1][mp.m_idx_y1:mp.m_idx_y2, mp.m_idx_x1:mp.m_idx_x2]))

x = []
for i in all_mesh:
    x.append(np.mean(i[mp.m_idx_y1:mp.m_idx_y2, mp.m_idx_x1:mp.m_idx_x2]))
y = [i for i in range(len(x))]

print(x[-1])

#plt.figure(2)
#plt.plot(y, x, "--", c='r')
#
#plt.figure(3)
#final_mesh = np.zeros((mp.rows, mp.cols))
#mp.update_nonboundaries(all_mesh[-1], final_mesh)
#
#masked = np.ma.masked_where(final_mesh < 0.01, final_mesh)
#ax = plt.subplot(111)
#im = ax.imshow(masked, cmap='rainbow')
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.xlabel('$T$ / K', labelpad=20)
#plt.colorbar(im, cax=cax)
