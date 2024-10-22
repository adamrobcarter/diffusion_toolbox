import common
import van_hove.van_hove as van_hove
import sys, time
import tqdm
import numpy as np
import matplotlib.pyplot as plt

def go(file, num_r_bins, max_r, outputfilename):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']

    num_timesteps = int(particles[:, 2].max()) + 1

    # num_r_bins = 100
    # max_r = 12
    num_time_origins = 20 # computation time linear in this

    gs = np.full((num_timesteps, num_r_bins), np.nan)

    
    # fig, ax = plt.subplots(1, 1)
    # ax.set_ylim(0, 2)

    n = num_timesteps # if file in ['eleanor0.01', 'alice0.02'] else 20

    time_origins = [int(i) for i in np.linspace(0, num_timesteps-1, num_time_origins)]

    for time_origin_index, time_origin in enumerate(tqdm.tqdm(time_origins)):
        # t0 = time.time()
        particles_at_t = particles[:, 2] == time_origin
        assert particles_at_t.sum() > 0
        # t1 = time.time()
        r_bin_edges, g, avg_density = van_hove.density_correlation(particles[particles_at_t, :], particles[particles_at_t, :], max_r=max_r, num_r_bins=num_r_bins, crop_border=10)
        gs[time_origin_index, :] = g
        # t2 = time.time()

        # a = t1-t0
        # b = t2-t1
        # t = t2-t0
        # print(f'{a/t:.2f}, {b/t:.2f}')

        # ax.clear()
        # ax.set_ylim(0, 3)

        # ax.errorbar(r_bin_edges, np.nanmean(gs, axis=0), yerr=np.nanstd(gs, axis=0)/np.sqrt(n), marker='.', linestyle='none')
        # ax.semilogy()

    common.save_data(outputfilename,
        g=gs, r=r_bin_edges,
        particle_diameter=data.get('particle_diameter'),
        num_time_origins=num_time_origins, num_r_bins=num_r_bins,
        pack_frac_given=data.get('pack_frac_given')
    )

    
for file in sys.argv[1:]:
    outputfilename = f'van_hove/data/g_{file}'
    num_r_bins = 300 # computation time (probably) linear in this
    max_r = 30 # computation time independent of this
    go(file, num_r_bins, max_r, outputfilename)
    
    # outputfilename = f'van_hove/data/g_{file}_a'
    # num_r_bins = 100 # computation time (probably) linear in this
    # max_r = 12 # computation time independent of this
    # go(file, num_r_bins, max_r, outputfilename)
    
    # outputfilename = f'van_hove/data/g_{file}_b'
    # num_r_bins = 200 # computation time (probably) linear in this
    # max_r = 12 # computation time independent of this
    # go(file, num_r_bins, max_r, outputfilename)