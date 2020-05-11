#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: MatteoRaso
"""

# Similar to the projectile function, this function simulates the trajectory
# of a projectile. Unlike the other projectile function, this function accounts
# for drag, which makes the results more realistic.

# Inputs:
# velocity of projectile relative to fluid as a 2D array
# density of fluid
# reference area of projectile
# drag coefficent
# initial position of the projectile as a 2D array
# gravity as a negative number

# Outputs:
# An array showing the projectile's trajectory throughout the simulation
# A visualization of the array

import numpy as np
import matplotlib.pyplot as plt


def drag_projectile(velocity, density, reference_area, drag_coefficent, position, gravity):
    if len(velocity) != 2 or len(position) != 2:
        print("Invalid input.")

    else:
        i = 0
        projection = np.array(position)

        while (position[1]) > 1e-4 or i < 1:
            # We need the angle for the drag force to act counter to the projectile.
            angle = np.tan(velocity[1] / velocity[0])
            drag_force = 0.5 * density * drag_coefficent * \
                (velocity[0] ** 2 + velocity[1] ** 2) * reference_area
            # Adding direction to the drag force.
            drag_force = np.array([np.cos(angle), np.sin(angle)]) * drag_force
            position[0] += velocity[0] * 0.01
            position[1] += velocity[1] * 0.01 + 0.5 * gravity * (0.01 ** 2)
            velocity[1] += gravity * 0.01
            velocity -= drag_force * 0.01
            projection = np.vstack((projection, position))
            i += 1

        plt.plot(projection[:, 0], projection[:, 1])
        return projection
