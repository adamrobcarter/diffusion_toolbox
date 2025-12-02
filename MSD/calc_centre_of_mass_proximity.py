import common
import MSD.MSD
import numpy as np
import tqdm

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particles     = data['particles']
        time_step     = data['time_step']
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']

        use_incremental = particles.size > 1e8
        use_incremental = True

        num_timesteps = int(particles[:, 2].max()) + 1
        av_particles_per_frame = particles.shape[0]/num_timesteps
        print('av particles per frame =', av_particles_per_frame)

        density = particles.shape[0] / (num_timesteps * particles[:, 0].max() * particles[:, 1].max())

        msds = []
        msd_uncs = []

        max_groupsize = av_particles_per_frame
        # if file == 'eleanor0.01': max_groupsize = 113
        # groupsizes = [int(N) for N in np.logspace(0, np.log10(np.sqrt(max_groupsize)), 40)**2]
        target_Ls = np.logspace(np.log10(2*2.82), np.log10(120*2.82), 50)
        box_sizes = [int(N) for N in target_Ls**2 * density]
        box_sizes = np.logspace(np.log10(10), np.log10(300), 8)

        msds     = np.full((len(box_sizes)), np.nan)
        msd_uncs = np.full_like(msds, np.nan)

        num_time_origins = 50

        # precompile
        # MSD.MSD.calc_centre_of_mass_onepoint_incremental_numba(particles, 1, num_time_origins)

        dframes = common.exponential_integers(1, num_timesteps-1, 100)

        # MSD.MSD.find_particles_at_frame(particles, num_timesteps)
        msds, num_used_particles = MSD.MSD.calc_centre_of_mass_proximity(particles, num_timesteps, window_size_x, window_size_y, num_time_origins, box_sizes)

        # msds, msd_uncs = MSD.MSD.calc_centre_of_mass_incremental_numba(particles, groupsizes, dframes, num_time_origins)

        common.save_data(f'MSD/data/msd_centre_of_mass_proximity_{file}',
                        msds=msds, msd_uncs=np.zeros_like(msds), box_sizes=box_sizes,
                        num_used_particles=num_used_particles,
                        time_step=data['time_step'],
                        density=density, num_time_origins=num_time_origins,
                        particle_diameter=data.get('particle_diameter'),
                        pixel_size=data['pixel_size'], window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y')
        )