import common
import numpy as np

if __name__ == "__main__":
    # test drift removal
    # num_timesteps = 1000
    # num_particles = 100
    # particles = []
    # drift_x = 0.03
    # drift_y = -0.05
    # for particle in range(num_particles):
    #     for timestep in range(num_timesteps):
    #         x = particle*10 + np.random.random() + drift_x*timestep
    #         y = particle*10 + np.random.random() + drift_y*timestep
    #         #   ^^^^^^^^^^^ this is so each particle has a different origin
    #         particles.append([x, y, timestep, particle])

    # particles = np.array(particles)
    # common.remove_drift(particles)

    # periodic unwrap
    particles = np.array([[0.5+i, 2, i, 0] for i in range(10)])
    print(particles)
    particles_wrapped = np.copy(particles)
    particles_wrapped[:, 0] = particles_wrapped[:, 0] % 4
    particles_unwrapped = common.periodic_unwrap(particles_wrapped, 2, [0, 1], [4, 4])
    print(particles_unwrapped)
    assert np.all(particles_unwrapped == particles)