#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:19:43 2022

@author: dallaire
"""

import numpy as np
from subfunctions import *
from define_rovers import *
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

rover, planet = define_rover_1()
Crr = 0.2
slope_list_deg = np.linspace(-10,35,25)
omega_max = np.zeros(len(slope_list_deg), dtype = float)
omega_nl = rover['wheel_assembly']['motor']['speed_noload']

# find where F_net == 0
for ii in range(len(slope_list_deg)):
    fun = lambda omega: F_net(omega, float(slope_list_deg[ii]), rover, planet, Crr)
    sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
    omega_max[ii] = sol.root
    
