import common
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        assert data['time_step'] == 1
        particles = data['particles']
        
        times = np.unique(particles[:, 2])
        long_timestep = times[1]
        print(f'keeping particles with time a multiple of {long_timestep}')
        to_keep = particles[:, 2] % long_timestep == 0
        print(f'keeping {to_keep.sum()/to_keep.size}')

        print('copying dictionary')
        newdata = dict(data)

        newdata['particles'][:, 2] /= long_timestep
        newdata['time_step'] = long_timestep
        del data
        
        times = np.unique(newdata['particles'][:, 2])
        long_timestep_new = times[1]
        assert long_timestep_new == 1, long_timestep_new

        filebase = file.split('_mixt')[0]
        common.save_data(f'particle_detection/data/particles_{filebase}_unmix.npz', **newdata)