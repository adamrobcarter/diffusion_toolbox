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

    n = num_timesteps if file in ['eleanor0.01', 'alice0.02'] else 20
    for timestep in tqdm.trange(n):
        particles_at_t = particles[:, 2] == timestep
        r_bin_edges, g, avg_density = van_hove.density_correlation(particles[particles_at_t, :], particles[particles_at_t, :], max_r=12, num_r_bins=num_r_bins, crop_border=10)
        gs[timestep, :] = g

        ax.clear()
        ax.set_ylim(0, 3)

        ax.errorbar(r_bin_edges, np.nanmean(gs, axis=0), yerr=np.nanstd(gs, axis=0)/np.sqrt(n), marker='.', linestyle='none')
        # ax.semilogy()

    common.save_data(f'van_hove/data/g_{file}', g=gs, r=r_bin_edges,
                     particle_diameter=data.get('particle_diameter'))

    