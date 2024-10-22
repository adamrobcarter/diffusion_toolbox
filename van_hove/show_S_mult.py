import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt

RESCALE_X_BY_DIAMETER = True
    
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

def go(file):
    data = common.load(f'van_hove/data/g_{file}.npz')
    g = data['g']
    r = data['r']
    particle_diameter = data.get('particle_diameter')
    
    g_mean = np.nanmean(g, axis=0)
    print('g', g_mean[:5])
    print('r', r.min(), r.max(), r.size)

    k, S = common.fourier(r, g_mean)

    safe_indexes = r/particle_diameter > 0.5
    r_safe = r     [safe_indexes]
    g_safe = g_mean[safe_indexes]
    k2, S2 = common.fourier(r_safe, g_safe)

    g_nodelta = g_mean[:]
    g_nodelta[0] = 0
    k3, S3 = common.fourier(r, g_nodelta)

    g4 = g_mean[:]
    g4[0] = 1
    k4, S4 = common.fourier(r, g4)
    # ax.set_ylabel('$g(r)$')
    # ax.set_ylim(0, 3)

    if RESCALE_X_BY_DIAMETER:
        x  = k  * particle_diameter
        x2 = k2 * particle_diameter
        x3 = k3 * particle_diameter
        x4 = k4 * particle_diameter
        ax.set_xlabel('$k\sigma$')
    else:
        x  = k
        x2 = k2
        x3 = k3
        x4 = k4
        ax.set_xlabel('$k$')
    
    # ax.tick_params(axis='y', colors='tab:blue')
    
    # ax.errorbar(x, S, marker='o', linestyle='none', color='tab:blue')
    
    # ax.errorbar(x3[1:], S3[1:], marker='o', linestyle='none', color='tab:cyan')
    
    ax.errorbar(x4[1:], S4[1:], marker='o', linestyle='-', label=file)
    

    # ax2 = ax.twinx()

    # ax2.errorbar(x2[1:], S2[1:], marker='o', linestyle='none', color='tab:orange')
    # ax2.tick_params(axis='y', colors='tab:orange')
    # ax.semilogy()


    # ax.errorbar(r/particle_diameter, g_mean, yerr=np.nanstd(g, axis=0)/np.sqrt(g.shape[0]), marker='o', color='tab:green')#, linestyle='none')
    
    # print(r_safe[max_index]/particle_diameter, g_safe[max_index])
    # ax.scatter(r_safe[max_index]/particle_diameter, g_safe[max_index], color='red', zorder=10)

    # ax.semilogy()

    

if __name__ == '__main__':
    for file in common.files_from_argv('van_hove/data', 'g_'):
        go(file)
    
    ax.legend()
    ax.hlines(1, *ax.get_xlim(), color='grey')
    common.save_fig(fig, f'van_hove/figures_png/S_of_k_mult.png')
