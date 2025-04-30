import common
import matplotlib.pyplot as plt
import numpy as np
import particle_linking.potential_at_z
import scipy

fig, ax = plt.subplots(1, 1)

for i, file in enumerate(files := common.files_from_argv('particle_linking/data/', 'trajs_')):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    a = data['particle_diameter']/2
    z = particles[:, 2]

    bins = np.linspace(1, 1.6, 160)
    bin_centers = (bins[1:] + bins[:-1])/2
    
    n_blobs = data['num_blobs']
    T = data['T']
    method = data['method']
    dt = data['dt']
    counts, bin_edges, _ = ax.hist(z/a, facecolor=common.tab_color(i), density=True, bins=bins, alpha=0.5, label=f'simulation, nblobs = {n_blobs}, T={T}K, {method}, dt={dt}')
    # counts, _ = np.histogram(z/a, bins=bins)
    # print(count)
    histogram_area = np.sum(counts * np.diff(bins))
    counts_normalised = counts/histogram_area
    # ax.bar(bin_centers, counts_normalised, width=np.diff(bins)[0], alpha=0.5, label=f'nblobs = {n_blobs}, T={T}K, {method}', color='red')
    # print('hist area', histogram_area)

    
    # z_th = np.linspace(1, 1.6, 1000) * data['particle_diameter']/2
    z_th = np.linspace(1, 1.6, 1000)
    z_over_a_th = z_th * data['particle_diameter']/2
    n_blobs = data['num_blobs']
    U = particle_linking.potential_at_z.potential_at_z(z_over_a_th, data['particle_diameter']/2, data['debye_length'], data['repulsion_strength'], n_blobs, data['T'], data['bare_particle_diameter']/2)

    atto = 1e-18
    k_B = scipy.constants.k / atto
    kT = k_B * data['T'] # in aJ
    probability = np.exp(-U / kT)
    assert probability.sum() > 0
    normalisation = scipy.integrate.trapezoid(probability, z_th)
    # normalisation = np.sum(probability * np.diff(z_th)[0])
    # normalisation = np.sum(probability * np.diff(bins))
    assert normalisation > 0
    probability_normalised = probability / normalisation
    ax.plot(z_over_a_th/a, probability_normalised, label=f'theory, nblobs = {n_blobs}', color=common.tab_color(i))
    
    # print('diff', np.mean(counts_normalised - probability_normalised))

    print('integral', scipy.integrate.trapezoid(probability, z_over_a_th))

    common.save_data('hist_theory.npz', z=z_over_a_th, P=probability)

ax.set_xlabel('$z/a$')
ax.legend()

filename = '_'.join(files)
common.save_fig(fig, f'particle_linking/figures_png/hist_z_{file}.png')