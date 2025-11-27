import common
import MSD.MSD
import numpy as np
import tqdm

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particles = data['particles']
        time_step = data['time_step']

        use_incremental = particles.size > 1e8
        use_incremental = True

        num_timesteps = particles[:, 2].max()
        av_particles_per_frame = particles.shape[0]/num_timesteps
        print('av particles per frame =', av_particles_per_frame)

        msds = []
        msd_uncs = []

        if not use_incremental:
            data_ = MSD.MSD.reshape(particles)
        groupsizes = [int(N) for N in np.logspace(0, np.log10(np.sqrt(av_particles_per_frame)), 10)**2]

        msds     = np.full((len(groupsizes)), np.nan)
        msd_uncs = np.full_like(msds, np.nan)

        density = particles.shape[0] / (num_timesteps * particles[:, 0].max() * particles[:, 1].max())
        print(f'density = {density:.2f}')
        print(4/np.pi * 0.02/2.8**2, 4/np.pi * 0.34/2.8**2)

        num_time_origins=25

        # precompile
        # MSD.MSD.calc_centre_of_mass_onepoint_incremental_numba(particles, 1, num_time_origins)

        # progress = tqdm.tqdm(groupsizes)
        if use_incremental:
            msds, msd_uncs = MSD.MSD.calc_centre_of_mass_onepoint_incremental(particles, groupsizes, num_time_origins)
        else:
            assert True, 'need to update nonincremental with new loop structure'
            msd, msd_unc = MSD.MSD.calc_centre_of_mass_onepoint(data_, groupsizes, num_time_origins)

        common.save_data(f'MSD/data/msd_centre_of_mass_{file}',
                        msds=msds, msd_uncs=msd_uncs, groupsizes=groupsizes,
                        time_step=data['time_step'], density=density,
                        num_time_origins=num_time_origins, particle_diameter=data.get('particle_diameter'),
                        pixel_size=data['pixel_size'], window_size_x=data['window_size_x'], window_size_y=data['window_size_y'])