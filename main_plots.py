#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 04:15:50 2018

@author: ricktjwong
"""

import matplotlib.pyplot as plt
import modules.plots.microprocessor_ceramic_plots as mcp
import modules.plots.natural_convection_sweep as ncs
import modules.plots.forced_convection_sweep as fcs
import modules.plots.simulated_annealing_analysis as sa

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
plt.rcParams['lines.linewidth'] = 2
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['figure.figsize'] = (15, 10)


"""
Question 3:
No Heat Sink
Ceramic-Microprocessor cooling
"""

mcp.plot_figures()

"""
Question 4:
Heat Sink with Natural Convection
"""

ncs.plot_sweeps()
sa.plot_action_costs()


"""
Question 5
Heat Sink with Forced Convection
Estimated total run time is 2 hours
"""

fcs.plot_sweeps()

plt.show()
