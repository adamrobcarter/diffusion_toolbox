import common
import numpy as np

for file in common.files_from_argv('preprocessing/data', 'particles'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    # particles = data['particles']
    # times, particles_at_time = np.unique(particles[:, 2], return_counts=True)
    # print(particles_at_time.max(), particles_at_time.mean())