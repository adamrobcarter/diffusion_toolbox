import common
import sys
import numpy as np
import matplotlib.pyplot as plt

#### this should reuse code from show_g



def go(file, ax, color=None):
    data = common.load(f'van_hove/data/g_{file}.npz')
    g = data['g']
    r = data['r']
    particle_diameter = data.get('particle_diameter')
    print('r max', r.max())
    
    ax.set_ylabel(r'$g(r)$')
    ax.set_xlabel(r'$r/a$')
    ax.set_ylim(0, 3)
    
    g_mean = np.nanmean(g, axis=0)
    print(g_mean[r/particle_diameter>0.5])

    a = particle_diameter/2

    safe_indexes = r/particle_diameter > 0.5
    r_safe = r     [safe_indexes]
    g_safe = g_mean[safe_indexes]
    max_index = np.nanargmax(g_safe)
    # print(g)
    print('max g', r_safe[max_index])

    if 'pack_frac' in data:
        label = fr'$\phi={data["pack_frac"]:.2f}$'
    else:
        label = file


    ax.errorbar(r/a, g_mean, yerr=np.nanstd(g, axis=0)/np.sqrt(g.shape[0]), marker='o', label=label,
                color=color)#, linestyle='none')
    
    # ax.scatter(r_safe[max_index]/a, g_safe[max_index], color='red', zorder=10)


if __name__ == '__main__':
    fig, ax = plt.subplots()

    for i, file in enumerate(sys.argv[1:]):
        color = plt.cm.copper(i/(len(sys.argv[1:]))-0.5)
        go(file, ax, color=color)

    ax.set_xlim(0, 10)

    # ax.semilogy()
    ax.legend()
    ax.axhline(1.0, color='gray', linestyle=':')
    common.save_fig(fig, f'van_hove/figures_png/g_of_r_{file}.png')