import common
import sys
import numpy as np
import matplotlib.pyplot as plt

def go(file, ax):
    data = common.load(f'van_hove/data/g_{file}.npz')
    g = data['g']
    r = data['r']
    particle_diameter = data.get('particle_diameter')
    print('r max', r.max())
    
    ax.set_ylabel('$g(r)$')
    ax.set_xlabel(r'$r/\sigma$')
    ax.set_ylim(0, 3)
    
    g_mean = np.nanmean(g, axis=0)
    print(g_mean[r/particle_diameter>0.5])

    safe_indexes = r/particle_diameter > 0.5
    r_safe = r     [safe_indexes]
    g_safe = g_mean[safe_indexes]
    max_index = np.nanargmax(g_safe)
    # print(g)
    print('max g', r_safe[max_index])

    ax.errorbar(r/particle_diameter, g_mean, yerr=np.nanstd(g, axis=0)/np.sqrt(g.shape[0]), marker='o', label=file)#, linestyle='none')
    
    print(r_safe[max_index]/particle_diameter, g_safe[max_index])
    ax.scatter(r_safe[max_index]/particle_diameter, g_safe[max_index], color='red', zorder=10)


if __name__ == '__main__':
    fig, ax = plt.subplots()

    for file in sys.argv[1:]:
        go(file, ax)

    # ax.semilogy()
    ax.legend(fontsize=6)
    ax.axhline(1.0, color='gray', linestyle=':')
    common.save_fig(fig, f'van_hove/figures_png/g_of_r_{file}.png')