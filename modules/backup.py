#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:57:54 2018

@author: ricktjwong
"""

    all_mesh = [l_mesh, r_mesh, u_mesh, d_mesh]
    for i in range(len(l_mesh[:,0])):
        if l_mesh[:,0][i] < 20:
            print(l_mesh[:,0][i])
            l_mesh[:,0][i] = 20
        else:
            l_mesh[i][0] = l_mesh[i][-1] - c * (l_mesh[i][1] - T_a) ** (4/3)
    for i in range(len(r_mesh[:,-1])):
        if r_mesh[:,-1][i] < 20:
            r_mesh[:,-1][i] = 20
        else:
            r_mesh[i][-1] = r_mesh[i][1] - c/2 * (r_mesh[i][1] - T_a) ** (4/3)
    for i in range(len(u_mesh[0])):
        if u_mesh[0][i] < 20:
            u_mesh[0][i] = 20
        else:
            u_mesh[:,i][0] = u_mesh[:,i][1] - c/2 * (u_mesh[:,i][1] - T_a) ** (4/3)
    for i in range(len(d_mesh[-1])):
        if d_mesh[-1][i] < 20:
            d_mesh[-1][i] = 20    
        else:
            d_mesh[:,i][-1] = d_mesh[:,i][1] - c/2 * (d_mesh[:,i][1] - T_a) ** (4/3)
    for k in all_mesh:
        for j in range(len(k[0])):
            for i in range(len(k)):
                if k[i][j] < 20:
                    k[i][j] = 20