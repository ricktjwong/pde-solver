#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:55:58 2018
@author: ricktjwong
"""

import numpy as np
from scipy import signal

k_m = 150                   # Conductivity of silicon Microchip in W/m K
k_c = 230                   # Conductivity of ceramic block in W/m K
k_a = 248                   # Conductivity of aluminium fin in W/m K
T_a = 20                    # Ambient temperature

class HeatStructure():
    def __init__(self, scale, b, c, f_h, n_fins,
                 conv_ratio = 1E-6, convection_type = "forced"):
        self.scale, self.b, self.c, self.f_h, self.n_fins, self.conv = \
        scale, b, c, f_h, n_fins, conv_ratio
        # Step size h (in m)
        self.h = 1E-3 / self.scale
        # Change in temperature of microprocessor every s
        self.rho = 0.25 * self.h ** 2 * 500 * 1E6 / k_m
        if convection_type == "natural":
            self.phi_m = self.h * 2.62 / k_m 
            self.phi_c = self.h * 2.62 / k_c
            self.phi_a = self.h * 2.62 / k_a
            self.power = 4 / 3
        elif convection_type == "forced":
            self.phi_m = self.h * 2 / k_m * (11.54 + 5.7 * 20)
            self.phi_c = self.h * 2 / k_c * (11.54 + 5.7 * 20)
            self.phi_a = self.h * 2 / k_a * (11.54 + 5.7 * 20)
            self.power = 1            
        self.T = self.b + self.c
        self.rows = (self.f_h + 7) * self.scale + 2
        self.cols = (self.T*(self.n_fins-1) + self.c) * self.scale + 2
        self.mesh = np.zeros((self.rows, self.cols)) 
        self.initialise_indices()
        self.initialise_temp()
        self.initialise_boundaries()
    
    def initialise_indices(self):
        """
        Initialise indices
        """        
        # Initialise row indices for fins
        self.f_idx_y1, self.f_idx_y2 = 1, (self.f_h*self.scale)+1
        # Initialise row and column indices for fin block
        self.fb_idx_x1, self.fb_idx_x2 = 1, self.cols - 1
        self.fb_idx_y1, self.fb_idx_y2 = \
        self.f_idx_y2, self.f_idx_y2 + 4 * self.scale
        # Initialise row and column indices for ceramic block
        self.c_idx_x1 = int((self.cols - 20*self.scale)/2)
        self.c_idx_x2 = self.c_idx_x1 + 20*self.scale
        self.c_idx_y1, self.c_idx_y2 = \
        self.fb_idx_y2, self.fb_idx_y2 + 2 * self.scale
        # Initialise row and column indices for microprocessor block
        self.m_idx_x1 = self.c_idx_x1 + 3 * self.scale
        self.m_idx_x2 = self.m_idx_x1 + 14 * self.scale
        self.m_idx_y1, self.m_idx_y2 = \
        self.c_idx_y2, self.c_idx_y2 + 1 * self.scale    

    def initialise_temp(self):
        """
        Set start temperature of the components to ambient temperature T_a
        """
        # Initialise fin mesh points to T_a
        for i in range(self.n_fins):
            self.mesh[self.f_idx_y1:self.f_idx_y2,
                      (self.T*i*self.scale+1):(self.c+self.T*i) \
                      * self.scale+1] = T_a
        # Initialise the fin block mesh points to T_a        
        self.mesh[self.fb_idx_y1:self.fb_idx_y2,
                  self.fb_idx_x1:self.fb_idx_x2] = T_a
        # Initialise the ceramic block mesh points to T_a
        self.mesh[self.c_idx_y1:self.c_idx_y2,
                  self.c_idx_x1:self.c_idx_x2] = T_a
        # Initialise the microprocessor mesh points to T_a
        self.mesh[self.m_idx_y1:self.m_idx_y2,
                  self.m_idx_x1:self.m_idx_x2] = T_a
        
    def initialise_boundaries(self):
        c_mesh, m_mesh, fb_mesh, f_mesh = self.update_all_boundaries(self.mesh)
        self.update_mesh(self.mesh, c_mesh, m_mesh, fb_mesh, f_mesh)
    
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
        l_mesh[:,0] = l_mesh[:,-1] - c * (l_mesh[:,1] - T_a) ** self.power
        r_mesh[:,-1] = r_mesh[:,0] - c * (r_mesh[:,1] - T_a) ** self.power
        u_mesh[0] = u_mesh[-1] - c * (u_mesh[1] - T_a) ** self.power
        d_mesh[-1] = d_mesh[0] - c * (d_mesh[1] - T_a) ** self.power
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
        fb_mesh = m[self.fb_idx_y1-1:self.fb_idx_y2+1,
                    self.fb_idx_x1-1:self.fb_idx_x2+1].copy()
        self.update_boundaries(fb_mesh, self.phi_a)
        f_mesh = [m[self.f_idx_y1-1:self.f_idx_y2+1,
                    (self.T*i*self.scale):(self.c+self.T*i)
                    * self.scale+2].copy() for i in range(self.n_fins)]
        for f in f_mesh:
            self.update_boundaries(f, self.phi_a)
        return c_mesh, m_mesh, fb_mesh, f_mesh
    
    def update_mesh(self, m, c_mesh, m_mesh, fb_mesh, f_mesh):
        """
        Update original mesh with the boundary values of different components. 
        The components are stacked in an order that preserves the correct
        boundaries. In this case we stack in the order: fin block, fins,
        ceramic, microchip
        """
        # Fin block with all boundaries
        m[self.fb_idx_y1-1:self.fb_idx_y2+1,
          self.fb_idx_x1-1:self.fb_idx_x2+1] = fb_mesh.copy()
        # Fins without the bottom boundary
        for i in range(self.n_fins):
            m[self.f_idx_y1-1:self.f_idx_y2,
              (self.T*i*self.scale):(self.c+self.T*i)*self.scale+2] \
              = f_mesh[i][0:-1].copy()
        # Handle the overlap boundaries
        m[self.fb_idx_y1-1][self.c*self.scale+1:-1][::self.T*self.scale] = \
        (m[self.fb_idx_y1-1][self.c*self.scale+1:-1][::self.T*self.scale] + \
         fb_mesh[0][self.c*self.scale+1:-1][::self.T*self.scale]) / 2
        m[self.fb_idx_y1-1][self.T*self.scale:][::self.T*self.scale] = \
        (m[self.fb_idx_y1-1][self.T*self.scale:][::self.T*self.scale] + \
         fb_mesh[0][self.T*self.scale:][::self.T*self.scale]) / 2
        # Ceramic without the top boundary
        m[self.c_idx_y1:self.c_idx_y2+1,
          self.c_idx_x1-1:self.c_idx_x2+1] = c_mesh[1:].copy()
        # Handle the overlap boundaries
        m[self.c_idx_y1][self.c_idx_x1 - 1] = \
        (m[self.c_idx_y1][self.c_idx_x1 - 1] + fb_mesh[-1][self.c_idx_x1 - 1]) / 2
        m[self.c_idx_y1][self.c_idx_x2] = \
        (m[self.c_idx_y1][self.c_idx_x2] + fb_mesh[-1][self.c_idx_x2]) / 2
#        # Microprocessor without the top boundary
        m[self.m_idx_y1:self.m_idx_y2+1,
          self.m_idx_x1-1:self.m_idx_x2+1] = m_mesh[2:].copy()
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
        # Update matrix for the fin block
        n[self.fb_idx_y1:self.fb_idx_y2, self.fb_idx_x1:self.fb_idx_x2] = \
        m[self.fb_idx_y1:self.fb_idx_y2, self.fb_idx_x1:self.fb_idx_x2]
        # Update matrix for the fins
        for i in range(self.n_fins):
            n[self.f_idx_y1:self.f_idx_y2, 
              (self.T*i*self.scale+1):(self.c+self.T*i)*self.scale+1] = \
            m[self.f_idx_y1:self.f_idx_y2,
              (self.T*i*self.scale+1):(self.c+self.T*i)*self.scale+1]
    
    def solve_mesh(self):
        """
        Iterate and solve for the temperature at every mesh point till the 
        change in average temperature of the microprocessor is below the
        convergence_ratio
        """
        n = 0
#        all_mesh = []
        while (True):
            update = self.mesh.copy()
            # Define Jacobi filter kernel for convolution                                         
            kernel = np.array([[0,1,0],
                               [1,0,1],
                               [0,1,0]]) 
            # Perform 2D convolution with input data and Jacobi filter kernel
            out = signal.convolve2d(self.mesh, kernel,
                                    boundary='wrap', mode='same')/kernel.sum()
            self.update_nonboundaries(out, update)
            update[self.m_idx_y1:self.m_idx_y2,
                   self.m_idx_x1:self.m_idx_x2] += self.rho
            m_mean_temp_1 = np.mean(self.mesh[self.m_idx_y1:self.m_idx_y2,
                                              self.m_idx_x1:self.m_idx_x2])
            norm_temp_1 = np.linalg.norm(self.mesh)
            norm_temp_2 = np.linalg.norm(update)
            if norm_temp_2 / norm_temp_1 - 1 < self.conv: break
            c_mesh, m_mesh, fb_mesh, f_mesh = \
            self.update_all_boundaries(update)
            update = self.update_mesh(update, c_mesh, m_mesh,
                                      fb_mesh, f_mesh).copy()
            self.mesh = update.copy()
            if n % 10000 == 0:
                print(n)
#            all_mesh.append(self.mesh)
            n += 1
#        return all_mesh
        return m_mean_temp_1, n
