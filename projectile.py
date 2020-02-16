#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 01:33:06 2020

@author: matteo
"""
#A simple projectile simulation

#Inputs:
#velocity as a 2D array
#gravity as a negative number
#position as a 2D array

#Outputs:
#An array showing the projectile's trajectory throughout the simulation
#A visualization of the array

import numpy as np
import matplotlib.pyplot as plt

def projectile(velocity, gravity, position):
    if len(velocity) != 2 or len(position) != 2:
        print("Invalid input.")
        
    else:
        
        i = 0
        projection = np.array(position)
        while (position[1]) > 1e-4 or i < 1:
            position[0] += velocity[0] * 0.01
            position[1] += velocity[1] * 0.01 + 0.5 * gravity * (0.01 ** 2)
            velocity[1] += gravity * 0.01
            projection = np.vstack((projection, position))
            i += 1
            
        plt.plot(projection[:, 0], projection[:, 1])
        return projection
    