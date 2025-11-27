import common
import numpy as np
import tqdm

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']

        print('x')
        common.term_hist(particles[:, 0])

        print('y')
        common.term_hist(particles[:, 1])

        print('particles per frame')
        common.term_hist(particles[:, 2])

        num_timesteps = particles[:, 2].max() + 1
        print('num timesteps:', num_timesteps)