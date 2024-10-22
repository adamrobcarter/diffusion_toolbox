import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt

def go(file):
    data = common.load(f'van_hove/data/g_{file}.npz')
    g = data['g']
    r = data['r']
    particle_diameter = data.get('particle_diameter')
    print('r max', r.max())
    
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    ax.set_ylabel('$g(r)$')
    ax.set_xlabel('$r/\sigma$')
    ax.set_ylim(0, 3)
    
    g_mean = np.nanmean(g, axis=0)
    print(g_mean[r/particle_diameter>0.5])

    safe_indexes = r/particle_diameter > 0.5
    r_safe = r     [safe_indexes]
    g_safe = g_mean[safe_indexes]
    max_index = np.nanargmax(g_safe)
    # print(g)
    print('max g', r_safe[max_index])

    ax.errorbar(r/particle_diameter, g_mean, yerr=np.nanstd(g, axis=0)/np.sqrt(g.shape[0]), marker='o', color='tab:green')#, linestyle='none')
    
    print(r_safe[max_index]/particle_diameter, g_safe[max_index])
    ax.scatter(r_safe[max_index]/particle_diameter, g_safe[max_index], color='red', zorder=10)

    # ax.semilogy()
    common.save_fig(fig, f'van_hove/figures_png/g_of_r_{file}.png')

    

if __name__ == '__main__':
    for file in sys.argv[1:]:
        go(file)