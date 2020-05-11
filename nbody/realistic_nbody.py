#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
@author: MatteoRaso
"""
# A more realistic n-body simulation.
# We will use it from the command-line.
# sys.argv[1] is a csv file containing the parameters of the simulation.
# The columns are time, box_width, box_height, box_depth
# Every other input is a csv file representing a particle.
# The columns are mass, x_position, y_position, z_position,  radius, x_velocity, y_velocity, z_velocity, charge
# The outputs are a unique csv file for each file.
# The columns are x_position, y_position, z_position, x_velocity, y_velocity, z_velocity, energy, momentum, current_time


import numpy as np
import sys
import csv
import auxillary as ax

epsilion = 1e-3
time_step = 0.01  # We will evaluate the simulation every 0.01 seconds.
G = 6.67e-11  # Gravitational constant
k = 8.99e9  # Electric constant


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
        current_particle.charge = float(row[8])
        f.close()
        current_particle.make_output_file()
        current_particle.init_array()
        particle_list.append(current_particle)

    current_time = 0
    updated_particle_list = particle_list[:]
    while time > current_time:
        current_time += time_step
        for particle in particle_list:
            ax.reflection(particle)
            # We remove the current particle so we don't accidentally have the particle collide with itself.
            updated_particle_list.remove(particle)
            for i in range(0, len(updated_particle_list)):
                second_particle = updated_particle_list[i]
                # The distance between the center of the particles
                distance = np.sqrt(sum((particle.position - second_particle.position))
                                   ** 2) + particle.radius + second_particle.radius

                if (np.sqrt(sum((particle.position - second_particle.position)) ** 2) - particle.radius - second_particle.radius) < epsilion:
                    ax.collision(particle, second_particle, current_time)
                    print("[*] Collision at " + str(current_time) + ".")

                # Updating the forces
                # We want each of the forces to be in the same direction as the second particle,
                # so we multiply the forces by the unit vector of the distance between the two.
                unit_vector = (second_particle.position - particle.position +
                               particle.radius + second_particle.radius) / distance
                particle.net_force += unit_vector * G * particle.mass * \
                    second_particle.mass / (distance ** 2)
                # We don't want to waste time calculating the electric force if it's 0.
                if particle.charge != 0 and second_particle.charge != 0:
                    # We want to use -= instead of += since the force having a positive value
                    # (i.e. sign(q_1) == sign(q_2) means that the force is repulsive.
                    # We could just multiply the value by -1 and keep the +=, but this
                    # is more efficent and elegant.
                    particle.net_force -= unit_vector.astype(
                        float) * k * particle.charge * second_particle.charge / (distance ** 2)

            particle.update_velocity(time_step)
            particle.update_energy()
            particle.update_momentum()
            current_array = np.array([particle.position[0], particle.position[1], particle.position[2], particle.velocity[0],
                                      particle.velocity[1], particle.velocity[2], particle.energy, particle.momentum, current_time])
            particle.position += particle.velocity * time_step
            particle.array = np.vstack((particle.array, current_array))

            updated_particle_list = particle_list[:]

    for particle in particle_list:
        particle.write_to_output_file()


if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) < 2:
    print("This is an n-body simulation that assumes no attraction or repulsion")
    print("between the particles. To use this program, give at least two csv")
    print("files giving the data for a particle. The data should be listed")
    print("as mass, velocity, radius, x_position, y_position, z_position.")

elif __name__ == '__main__':
    main(sys.argv[1:])
