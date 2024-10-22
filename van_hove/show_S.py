import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt
import countoscope_theory.structure_factor

RESCALE_X_BY_DIAMETER = True

def go(file):
    data = common.load(f'van_hove/data/g_{file}.npz')
    g = data['g']
    r = data['r']
    particle_diameter = data.get('particle_diameter')
    
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    
    g_mean = np.nanmean(g, axis=0)
    print('g', g_mean[:5])
    print('r step', r[1]-r[0])
    print('g[0]', g_mean[0])

    k, S = common.fourier(r, g_mean)

    safe_indexes = r/particle_diameter > 0.5
    r_safe = r     [safe_indexes]
    g_safe = g_mean[safe_indexes]
    k2, S2 = common.fourier(r_safe, g_safe)

    g_nodelta = g_mean[:]
    g_nodelta[0] = 0
    k3, S3 = common.fourier(r, g_nodelta)

    # g_mean = np.concatenate([g_mean, np.ones(g_mean.size*2)])

    g4 = g_mean[:]
    g4[0] = 1
    k4, S4 = common.fourier(r, g4)
    # ax.set_ylabel('$g(r)$')
    # ax.set_ylim(0, 3)
    
    g5 = g_mean[:] - 1
    g5[0] = 1000
    # g5[0] = 1
    k5, S5 = common.fourier(r, g5)
    S5 /= 1000
    k5 *= 2 * np.pi
    
    g6 = g_mean[:] - 1
    g6[0] = 100
    # g5[0] = 1
    k6, S6 = common.fourier(r, g6)
    S6 /= 100
    k6 *= 2 * np.pi
    
    g7 = g_mean[:] - 1
    a = 37
    g7[0] = a
    k7, S7 = common.fourier(r, g7)
    S7 /= a
    k7 *= 2 * np.pi
    
    k_th = np.linspace(0.5, 10, 100)

    if RESCALE_X_BY_DIAMETER:
        # x  = k  * particle_diameter
        # x2 = k2 * particle_diameter
        # x3 = k3 * particle_diameter
        # x4 = k4 * particle_diameter
        x5 = k5 * particle_diameter
        x6 = k6 * particle_diameter
        x7 = k7 * particle_diameter
        x_th = k_th * particle_diameter
        ax.set_xlabel('$k\sigma$')
    else:
        # x  = k
        # x2 = k2
        # x3 = k3
        # x4 = k4
        x5 = k5
        x6 = k6
        x7 = k7
        x_th = k_th
        ax.set_xlabel('$k$')
    
    
    # ax.errorbar(x, S, marker='o', linestyle='none', color='tab:blue')
    
    # ax.errorbar(x3[1:], S3[1:]+1, marker='none', linestyle='-', label='g[0]=0')
    
    # ax.errorbar(x4[1:], S4[1:]+1, marker='none', linestyle='-', label='g[0]=1')
    
    ax.plot(x_th, countoscope_theory.structure_factor.hard_spheres_2d(k_th, 0.34, 3.03), color='grey', label='$\sigma=3.03$')

    # # ax.tick_params(axis='y', colors='tab:blue')
    # # ax2 = ax.twinx()
    # # ax2.tick_params(axis='y', colors='tab:orange')

    # ax.errorbar(x2[1:], S2[1:]+1, marker='none', linestyle='-', label='g[0-0.5]=0')

    ax.errorbar(x5[1:], S5[1:], marker='none', linestyle='-', label='new 1000')
    ax.errorbar(x6[1:], S6[1:], marker='none', linestyle='-', label='new 100')
    ax.errorbar(x7[1:], S7[1:], marker='none', linestyle='-', label='new 10')
    
    # ax.semilogy()
    
    print('density', common.pack_frac_to_density(0.34, particle_diameter), 1/common.pack_frac_to_density(0.34, particle_diameter))

    # ax.errorbar(r/particle_diameter, g_mean, yerr=np.nanstd(g, axis=0)/np.sqrt(g.shape[0]), marker='o', color='tab:green')#, linestyle='none')
    
    # print(r_safe[max_index]/particle_diameter, g_safe[max_index])
    # ax.scatter(r_safe[max_index]/particle_diameter, g_safe[max_index], color='red', zorder=10)

    # ax.semilogy()
    ax.legend()
    ax.hlines(1, *ax.get_xlim(), color='grey')
    common.save_fig(fig, f'van_hove/figures_png/S_of_k_{file}.png', dpi=200)

    

if __name__ == '__main__':
    for file in common.files_from_argv('van_hove/data', 'g_'):
        go(file)