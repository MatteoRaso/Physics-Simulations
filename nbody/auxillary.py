#!usr/bin/python
#Contains all of the auxillary functions that we'll be needing in nbody.py

import numpy as np
import sys
import csv
import os
import pandas as pd

epsilion = 1e-3
norm = np.linalg.norm #Makes things more readable.

with open(sys.argv[1], 'r') as csvfile:
       reader = csv.reader(csvfile, delimiter = ',')
       print(reader)
       for row in reader:
           time = float(row[0])
           #Specifies the boundaries of the box where this simulation will be taking place.
           max_box = np.array([float(row[1]), float(row[2]), float(row[3])]) / 2
           min_box = -max_box

class particle():
    def init(self):
        self.mass = 0 #Kg
        self.velocity = np.array([0,0,0]) #m/s
        self.momentum = 0
        self.energy = 0
        self.radius = 0
        self.volume = 0
        self.position = np.array([0,0,0]) #Measured from centre.
        self.distance_from_centre = 0
        self.position_array = np.array([0,0,0])
        self.input_file = ''
        self.output_file = ''
        self.array = []
        self.charge = 0

    def update_momentum(self):
        self.momentum = self.mass * np.sqrt(sum(self.velocity ** 2)) #Pythagorean Theorem

    def update_energy(self):
        self.energy = 0.5 * self.mass * sum(self.velocity ** 2)

    def update_distance(self):
          self.distance_from_centre = np.sqrt(sum(self.position ** 2))

    def make_output_file(self):
        self.output_file = self.input_file.strip('.csv') + '_output' + '.csv'

    def write_to_output_file(self):
        header = ['x_position', 'y_position', 'z_position', 'x_velocity', 'y_velocity', 'z_velocity',
                       'energy', 'momentum', 'time']
        self.array = np.vstack((header, self.array))
        #We need every element of the array to be a string for savetxt to work properly.
        self.array = self.array.astype(str) 
        np.savetxt(self.output_file, self.array, delimiter=',', fmt = '%s')
    def init_array(self):
        self.array = np.array([self.position[0], self.position[1], self.position[2], self.velocity[0], self.velocity[1], self.velocity[2], self.energy, self.momentum, 0])

def collision(particle1, particle2, time):
    global num_collisions
    position1 = particle1.position
    position2 = particle2.position
    velocity1 = particle1.velocity
    velocity2 = particle2.velocity
    mass1 = particle1.mass
    mass2 = particle2.mass
    total_energy = particle1.energy + particle2.energy
    #We want to find the velocity component pointing directly to the other particle.
    position_vector_1 = position1 - position2
    position_vector_2 = -position_vector_1
    #For a quick understanding of the math, refer to Dan Flath's answer at
    #https://www.quora.com/How-do-I-find-the-vector-component-of-a-vector-in-the-direction-of-another-vector
    #In 3D, the equation (V -cD) * D = 0 expands to
    #V_xD_x - cD_x^2 +V_yD_y - cD_y^2 + V_zD_z - cD_z^2 = 0.
    #We are using the velocity as V and the position vector as D.
    #Since we know those values, we just have to rearrange the equation
    #to get our scaling factor.
    c1 = np.dot(velocity1, position_vector_1) / (sum(position_vector_1 ** 2))
    c2 = np.dot(velocity2, position_vector_2) / (sum(position_vector_2 ** 2))
    v_centre_1 = c1 * position_vector_1
    v_centre_2 = c2 * position_vector_2
    v_normal_1 = velocity1 - v_centre_1
    v_normal_2 = velocity2 - v_centre_2
    #We can now treat this as a 1D elastic collision.
    temp1 = v_centre_1
    temp2 = v_centre_2
    v_centre_1 = ((mass1 - mass2) / (mass1 + mass2)) * temp1 + 2 * (mass2 / (mass1 + mass2)) * temp2
    v_centre_2 = 2 * (mass1 / ( mass1 + mass2)) * temp1 + ((mass2 - mass1) / (mass1 + mass2)) * temp2
    velocity1 = v_normal_1 + v_centre_1
    velocity2 = v_normal_2 + v_centre_2
    particle1.velocity = velocity1
    particle2.velocity = velocity2
    particle1.update_momentum()
    particle1.update_energy()
    particle2.update_momentum()
    particle2.update_energy()
    if ((particle1.energy + particle2.energy - total_energy) < epsilion) != True:
          print("[!] Energy just increased at time " + str(time) + ". This defies physics!!!")

def reflection(particle):
    #Deals with a particle that gets too close to the edge of the box.
    for i in range(0, 3):
        if abs(particle.position[i] - max_box[i]) < epsilion or abs(particle.position[i] - min_box[i]) < epsilion:
             #Reverses the velocity component parallel to the side of the box the particle is about to hit.
             particle.velocity[i] *= -1

