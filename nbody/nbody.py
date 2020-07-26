#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
@author: MatteoRaso
"""
# A simple n-body simulation. It is assumed that there are no attractive
# or repulsive forces between the particles.
# We will use it from the command-line
# sys.argv[1] is a csv file containing the parameters of the simulation.
# The columns are time, box_width, box_height, box_depth
# Every other input is a csv file representing a particle.
# The columns are mass, x_position, y_position, z_position,  radius, x_velocity, y_velocity, z_velocity
# The outputs are a unique csv file for each file.
# The columns are x_position, y_position, z_position, x_velocity, y_velocity, z_velocity, energy, momentum, current_time
#
# Algorithm
#
# 1. read in the files
# 2. set the parameters
# 3. using a for-loop, initalize a particle object for every file other than the parameters file
# 4. in a while-loop, run reflection(particle) for every particle
# 4.1. check if each particle is close to any other particle
# 4.1.1. if a particle is close to another particle, run collision(particle1, particle2)
# 4.2. Update particle positions
# 4.2. save the particle properties and current_time to numpy array current_array
#4.3. particle.array = np.vstack((particle.array, current_array))
# 4.4. add time_step to current_time
# 4.5. go to 4.1. and repeat until current_time >= time
# 5. Use a for-loop and execute particle.write_to_output_file() for each particle

import numpy as np
import sys
import csv
import auxillary as ax

epsilion = 1e-3
time_step = 0.01  # We will evaluate the simulation every 0.01 seconds.


def main(*args):
    with open(sys.argv[1], 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            time = float(row[0])
            # Specifies the boundaries of the box where this simulation will be taking place.
            max_box = np.array(
                [float(row[1]), float(row[2]), float(row[3])]) / 2
            min_box = -max_box

    particle_list = []
    # Extracting the elements from the list
    # We don't want the first file since that's our parameters, not the particle attributes
    args = args[0][1:]
    for file in args:
        current_particle = ax.particle()
        current_particle.init()
        i = 0
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                current_particle.mass = float(row[0])
                current_particle.position = np.array(
                    [float(row[1]), float(row[2]), float(row[3])])
                current_particle.radius = float(row[4])
                current_particle.velocity = np.array(
                    [float(row[5]), float(row[6]), float(row[7])])

        current_particle.update_energy()
        current_particle.update_momentum()
        current_particle.update_distance()
        f = open(file, 'r')
        current_particle.input_file = f.name
        f.close()
        current_particle.make_output_file()
        current_particle.init_array()
        particle_list.append(current_particle)

    current_time = 0
    updated_particle_list = particle_list.copy()
    while time > current_time:
        current_time += time_step
        for particle in particle_list:
            ax.reflection(particle)
            # We remove the current particle so we don't accidentally have the particle collide with itself.
            updated_particle_list.remove(particle)
            for i in range(0, len(updated_particle_list)):
                second_particle = updated_particle_list[i]
                # The distance between the particles minus the radius of the particles
                if (np.sqrt(sum((particle.position - second_particle.position) ** 2)) - particle.radius - second_particle.radius) < epsilion:
                    ax.collision(particle, second_particle, current_time)
                    print("[*] Collision at " + str(current_time) + ".")

            current_array = np.array([particle.position[0], particle.position[1], particle.position[2], particle.velocity[0],
                                      particle.velocity[1], particle.velocity[2], particle.energy, particle.momentum, current_time])
            particle.position += particle.velocity * time_step
            particle.array = np.vstack((particle.array, current_array))

            updated_particle_list = particle_list.copy()

    for particle in particle_list:
        particle.write_to_output_file()


if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) < 2:
    print("This is an n-body simulation that assumes no attraction or repulsion")
    print("between the particles. To use this program, give at least two csv")
    print("files giving the data for a particle. The data should be listed")
    print("as mass, velocity, radius, x_position, y_position, z_position.")

elif __name__ == '__main__':
    main(sys.argv[1:])
