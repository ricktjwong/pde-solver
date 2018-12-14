#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 17:10:49 2018

@author: ricktjwong
"""

import numpy as np
import random
import utils.solvers as solv
import modules.heat_structure as hst
import time

"""
scale, b, c, f_h, n_fins
We vary b, c, f_h, n_fins
1 < b < 5
1 < c < 5
30 < f_h < 50
6 < n_fins < 50
"""

b_min, b_max = 1, 5
c_min, c_max = 1, 5
f_h_min, f_h_max = 30, 50
n_fins_min, n_fins_max = 6, 50


def get_action(idx):
    """
    Impose boundaries for min and max for the 4 different parameters
    """
    if idx == 0:
        if x0[0] == b_min:
            move = 1
        elif x0[0] == b_max:
            move = -1
        else:
            move = random.choice([-1, 1])
    elif idx == 1:
        if x0[1] == c_min:
            move = 1
        elif x0[1] == c_max:
            move = -1
        else:
            move = random.choice([-1, 1])            
    elif idx == 2:
        if x0[2] == f_h_min:
            move = 1
        elif x0[2] == f_h_max:
            move = -1
        else:
            move = random.choice([-1, 1])            
    elif idx == 3:
        if x0[3] == n_fins_min:
            move = 1
        elif x0[3] == n_fins_max:
            move = -1         
        else:
            move = random.choice([-1, 1])
    return move


def acceptance_probability(old, new, T):
    """
    If the new cost value is smaller than the old one, accept the new cost
    value with 100% probability. If the new cost value is large than the old,
    accept the new cost value with a probability equal to exp[-(new - old) / T]
    T gets smaller at every iteration, i.e. we accept bad cost values with
    lower probability
    """
    if new < old:
        a = 1
    else:
        a = np.exp((old - new) / T)
    return a

x0 = [2, 2, 30, 30]
min_cost = [1000, x0.copy()]
T = 1.0
T_min = 0.01
alpha = 0.9
costs = []
min_costs = []
x = []

while T > T_min:
    count = 0
    while(count < 100):
        b, c, f_h, n_fins = x0
        hs = hst.HeatStructure(2, b, c, f_h, n_fins, conv_ratio=1E-6,
                               convection_type="natural",
                               solver=solv.jacobi_solver)
        cost_new, n = hs.solve_mesh()
        ep = acceptance_probability(min_cost[0], cost_new, T)
        costs.append(cost_new)
        if ep > random.random():
            min_cost = cost_new, x0.copy()
            min_costs.append(min_cost)
        idx = random.choice([0, 1, 2, 3])
        move = get_action(idx)
        x0[idx] += move
        x.append(x0.copy())
        count += 1
    np.savetxt("temp_run" + str(T) + ".txt",
               [min_cost, min_costs, costs, x], fmt='%s')
    T = T * alpha
    
#x0 = [5, 5, 50, 50]
#b, c, f_h, n_fins = x0
#hs = hst.HeatStructure(2, b, c, f_h, n_fins, conv_ratio=1E-6,
#                       convection_type="natural", solver=solv.jacobi_solver)
#cost_new, n = hs.solve_mesh()
#print(cost_new)

