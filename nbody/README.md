A simple n-body simulation. It is assumed that there are no attractive or repulsive forces between the particles.
The first argument is a csv file containing the parameters of the simulation
The columns of that file are time, box_width, box_height, and box_depth
Every other input is a csv file representing a particle.
The columns of the particle csv files are mass, x_position, y_position, z_position,  radius, x_velocity, y_velocity, and z_velocity
The outputs are a unique csv file for each file.
The columns of the output csv files are  x_position, y_position, z_position, x_velocity, y_velocity, z_velocity, energy, momentum, and current_time
