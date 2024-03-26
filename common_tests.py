import common
import numpy as np

# test drift removal
num_timesteps = 1000
num_particles = 100
particles = []
drift_x = 0.03
drift_y = -0.05
for particle in range(num_particles):
    for timestep in range(num_timesteps):
        x = particle*10 + np.random.random() + drift_x*timestep
        y = particle*10 + np.random.random() + drift_y*timestep
        #   ^^^^^^^^^^^ this is so each particle has a different origin
        particles.append([x, y, timestep, particle])

particles = np.array(particles)
common.remove_drift(particles)