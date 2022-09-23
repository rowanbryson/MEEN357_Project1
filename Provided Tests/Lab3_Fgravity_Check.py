#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 15:07:01 2022

@author: dallaire
"""

from subfunctions import *
from define_rovers import *
import numpy as np

rover, planet = define_rover_1()
terrain_angle = np.array([-5.0, 0.0, 5.0, 10.0, 20.0, 30.0])  # degrees!
Fgt = F_gravity(terrain_angle, rover, planet)
print('Omega     Fgt')
for i in range(len(Fgt)):
    print('{:3.4F}    {:3.4F}'.format(terrain_angle[i], Fgt[i]))