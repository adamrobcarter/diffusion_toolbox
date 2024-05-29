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


    G = np.full((num_timesteps, num_r_bins), np.nan)

    progress = tqdm.tqdm(total=len(delta_ts)*n)

    for delta_t_index in range(len(delta_ts)):

        gs = np.full((num_timesteps, num_r_bins), np.nan)

        for timestep in range(n):
            particles_at_t0 = particles[:, 2] == timestep
            particles_at_t1 = particles[:, 2] == (timestep + delta_ts[delta_t_index])
            r_bin_edges, g, avg_density = van_hove.density_correlation(particles[particles_at_t0, :], particles[particles_at_t1, :], max_r=12, num_r_bins=num_r_bins, crop_border=10)
            gs[timestep, :] = g

            progress.update()

        G[delta_t_index, :] = np.nanmean(gs, axis=0)

    common.save_data(f'van_hove/data/G_{file}', G=G, r=r_bin_edges, t=delta_ts,
                     particle_diameter=data.get('particle_diameter'))

    