import common
import MSD.MSD
import numpy as np
import tqdm

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data/', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        time_step = data['time_step']

        use_incremental = particles.size > 1e8
        use_incremental = True

        num_timesteps = np.unique(particles[:, 2]).size
        av_particles_per_frame = particles.shape[0]/num_timesteps

        particles_per_frame = np.bincount(particles[:, 2].astype('int'))
        assert np.all(particles_per_frame == particles_per_frame[0])
        num_particles_per_frame = int(particles_per_frame[0])
        
        # need the data sorted by time
        print('sorting')
        particles = particles[particles[:, 2].argsort()]

        centre_of_mass = np.zeros((num_timesteps, 2))
        
        for t in tqdm.trange(num_timesteps):
            # t=int(t)
            particles_this_timestep = particles[t*num_particles_per_frame:(t+1)*num_particles_per_frame, :]
            assert np.all(particles_this_timestep[:, 2] == t)
            centre_of_mass[t, 0] = np.mean(particles_this_timestep[:, 0])
            centre_of_mass[t, 1] = np.mean(particles_this_timestep[:, 1])

        x_msd = MSD.MSD.msd_fft1d(centre_of_mass[:, 0])
        y_msd = MSD.MSD.msd_fft1d(centre_of_mass[:, 1])
        msd = x_msd + y_msd
        msd_unc = np.zeros_like(msd)

        # msd_sum    = np.zeros((num_timesteps))
        # n          = np.zeros((num_timesteps))
        # msd_sq_sum = np.zeros((num_timesteps)) 
        # for t_origin in tqdm.trange(num_timesteps):
        # # for t_origin in range(1000):
        #     msd_for_timeorigin = np.sum((centre_of_mass[t_origin:, :] - centre_of_mass[t_origin, :])**2, axis=1)
        #     msd_sum   [0:msd_for_timeorigin.size] += msd_for_timeorigin
        #     msd_sq_sum[0:msd_for_timeorigin.size] += msd_for_timeorigin**2
        #     n         [0:msd_for_timeorigin.size] += 1
        # msd = msd_sum / n
        # msd_unc = np.sqrt(msd_sq_sum / n - (msd_sum / n)**2)
        

        common.save_data(f'MSD/data/msd_centre_of_mass_entire_{file}',
                        msd=msd, msd_unc=msd_unc,
                        t = list(range(0, num_timesteps)) * data['time_step'],
                        #  density=density, num_time_origins=num_time_origins,
                        density=common.calc_density(particles, data['window_size_x'], data['window_size_y'], dimension=data.get('dimension', 2)),
                        N_particles=num_particles_per_frame,
                        particle_diameter=data.get('particle_diameter'),
                        pixel_size=data.get('pixel_size'),
                        window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y')
        )