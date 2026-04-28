import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt

for file in sys.argv[1:]:
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']

    num_timesteps = int(particles[:, 2].max())
    num_r_bins = 100

    gs = np.full((num_timesteps, num_r_bins), np.nan)

    
    fig, ax = plt.subplots(1, 1)
    ax.set_ylim(0, 2)
    
    delta_ts = range(0, 400, 10)

    n = num_timesteps - max(delta_ts) if file in ['eleanor0.01', 'alice0.02'] else 20


    G      = np.full((num_timesteps, num_r_bins), np.nan)
    G_self = np.full((num_timesteps, num_r_bins), np.nan)

    progress = tqdm.tqdm(total=len(delta_ts)*n)

    for delta_t_index in range(len(delta_ts)):

        gs = np.full((num_timesteps, num_r_bins), np.nan)
        gs_self = np.full((num_timesteps, num_r_bins), np.nan)

        for timestep in range(n):
            particles_at_t0 = particles[:, 2] == timestep
            particles_at_t1 = particles[:, 2] == (timestep + delta_ts[delta_t_index])
            r_bin_edges, g, avg_density = van_hove.density_correlation(particles[particles_at_t0, :], particles[particles_at_t1, :], max_r=12, num_r_bins=num_r_bins, crop_border=10)
            gs[timestep, :] = g

            # to do the self we need to know the identities of the particles
            common_ids = np.unique(np.concatenate((particles[particles_at_t0, 3], particles[particles_at_t1, 3])))
            # print(common_ids)
            # print(np.isin(particles[particles_at_t0, 3], common_ids).nonzero())
            # print(np.isin(particles[particles_at_t1, 3], common_ids).nonzero())
            common_particles_t0 = particles[particles_at_t0, :][np.isin(particles[particles_at_t0, 3], common_ids).nonzero()[0], :]
            common_particles_t1 = particles[particles_at_t1, :][np.isin(particles[particles_at_t1, 3], common_ids).nonzero()[0], :]
            print(common_particles_t0.shape, common_particles_t1.shape)
            assert common_particles_t0.shape == common_particles_t1.shape
            for i in range(common_particles_t0.shape[0]):
                assert common_particles_t0[i, 3] == common_particles_t1[i, 3]

            _, g_s, _ = van_hove.density_correlation_self(common_particles_t0, common_particles_t1, max_r=12, num_r_bins=num_r_bins, crop_border=10)
            gs_self[timestep, :] = g_s

            progress.update()

        G     [delta_t_index, :] = np.nanmean(gs, axis=0)
        G_self[delta_t_index, :] = np.nanmean(gs_self, axis=0)

    common.save_data(f'van_hove/data/G_{file}', G=G, G_self=G_self, r=r_bin_edges, t=delta_ts,
                     particle_diameter=data.get('particle_diameter'))

    