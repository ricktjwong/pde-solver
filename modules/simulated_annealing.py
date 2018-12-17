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

"""
scale, b, c, f_h, n_fins
b - fin separation
c - fin thickness
1 < b < 10
1 < c < 10
30 < f_h < 60
6 < n_fins < 60
"""

b_min, b_max = 1, 10
c_min, c_max = 1, 10
f_h_min, f_h_max = 30, 60
n_fins_min, n_fins_max = 6, 60


def get_action(idx):
    """
    Impose boundaries for min and max for the 4 different parameters
    """
    if idx == 0:
        # If we reach the minimum of parameter b, do a random choice of staying
        # put or taking one step forward
        if x0[0] == b_min:
            move = random.choice([0, 1])
        # If we reach the maximum of parameter b, do a random choice of staying
        # put or taking one step back            
        elif x0[0] == b_max:
            move = random.choice([0, -1])
        else:
            move = random.choice([-1, 1])
    elif idx == 1:
        if x0[1] == c_min:
            move = random.choice([0, 1])
        elif x0[1] == c_max:
            move = random.choice([0, -1])
        else:
            move = random.choice([-1, 1])            
    elif idx == 2:
        if x0[2] == f_h_min:
            move = random.choice([0, 1])
        elif x0[2] == f_h_max:
            move = random.choice([0, -1])
        else:
            move = random.choice([-1, 1])            
    elif idx == 3:
        if x0[3] == n_fins_min:
            move = random.choice([0, 1])
        elif x0[3] == n_fins_max:
            move = random.choice([0, -1])
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


def simulate_annealing(x0, T, T_min, alpha):
    min_costs = []
    min_actions = []
    min_cost = 1000
    while T > T_min:
        count = 0
        while(count < 100):
            b, c, f_h, n_fins = x0
            hs = hst.HeatStructure(2, b, c, f_h, n_fins, conv_ratio=1E-6,
                                   convection_type="natural",
                                   solver=solv.red_black_SOR)
            cost_new, n = hs.solve_mesh()
            ep = acceptance_probability(min_cost, cost_new, T)
            if ep > random.random():
                min_cost = cost_new
                min_action = x0.copy()
                min_costs.append(min_cost)
                min_actions.append(min_action)
            idx = random.choice([0, 1, 2, 3])
            move = get_action(idx)
            x0[idx] += move
            count += 1
        np.save("sim_annealing/costs" + str(T), min_costs)
        np.save("sim_annealing/action" + str(T), min_actions)
        T = T * alpha
        
x0 = [2, 2, 30, 30]
T = 1.0
T_min = 0.00001
alpha = 0.8

simulate_annealing(x0, T, T_min, alpha)
