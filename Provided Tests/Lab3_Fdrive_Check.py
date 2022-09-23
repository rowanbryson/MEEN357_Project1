#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 08:39:38 2022

@author: dallaire
"""

from subfunctions import *
from define_rovers import *
import numpy as np

rover, planet = define_rover_1()
omega = np.array([0.00, 0.50, 1.00, 2.00, 3.00, 3.80])
Fd = F_drive(omega, rover)
print('Omega      Fd')
for i in range(len(Fd)):
    print('{:3.4F}    {:3.4F}'.format(omega[i], Fd[i]))
    
    