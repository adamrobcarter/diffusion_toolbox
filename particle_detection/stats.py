# this file is just for seeing what's saved inside a data file

import common
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        time_column = data.get('dimension', 2)
        particles = data['particles']
        print(f'mean num particles, {particles.shape[0]/(particles[:, time_column].max()+1):.0f}')
        print(f'{particles[:, 0].min()} <= x <= {particles[:, 0].max()}')
        print(f'{particles[:, 1].min()} <= y <= {particles[:, 1].max()}')
        if data.get('dimension', 2) == 3:
            print(f'{particles[:, 2].min()} <= z <= {particles[:, 2].max()}')
            print(f'<z> = {particles[:, 2].mean()}')
        num_particles_at_timestep = np.bincount(particles[:, time_column].astype('int'))
        print('max particles at timestep', num_particles_at_timestep.max())
        print('min particles at timestep', num_particles_at_timestep.min())
        max_time_hours = data['max_time_hours']
        print(f't_max (from hours) = {max_time_hours*60*60:.2e}s = {max_time_hours*60*60:.0f}s')
        print(f'num timesteps', np.unique(particles[:, time_column]).size)