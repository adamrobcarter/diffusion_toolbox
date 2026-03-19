# this file is just for seeing what's saved inside a data file

import common
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']

        time_column = common.get_particles_column('t', data)
        x_column = common.get_particles_column('x', data)
        y_column = common.get_particles_column('y', data)
        print(f'mean num particles, {particles.shape[0]/(particles[:, time_column].max()+1):.0f}')
        print(f'{particles[:, x_column].min()} <= x <= {particles[:, x_column].max()}')
        print(f'{particles[:, y_column].min()} <= y <= {particles[:, y_column].max()}')
        if common.particles_has_z(data):
            z_column = common.get_particles_column('z', data)
            print(f'{particles[:, z_column].min()} <= z <= {particles[:, z_column].max()}')
            print(f'<z> = {particles[:, z_column].mean()}')

        num_particles_at_timestep = np.bincount(particles[:, time_column].astype('int'))
        print('max particles at timestep', num_particles_at_timestep.max())
        print('min particles at timestep', num_particles_at_timestep.min())

        max_time_hours = data['max_time_hours']
        print(f't_max (from hours) = {max_time_hours*60*60:.2e}s = {max_time_hours*60*60:.0f}s')
        print(f'num timesteps', np.unique(particles[:, time_column]).size)