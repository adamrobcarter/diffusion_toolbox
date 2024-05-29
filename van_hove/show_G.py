import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt

for file in sys.argv[1:]:
    data = common.load(f'van_hove/data/G_{file}.npz')
    G = data['G']
    r = data['r']
    t = data['t']
    particle_diameter = data.get('particle_diameter')
    
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    # ax.semilogy()

    def func(t_index):
        ax.clear()
        ax.set_ylabel('$G(r, \Delta t)$')
        ax.set_xlabel('$r/\sigma$')
        # ax.semilogy() set_ylim(0.5, 30)
        ax.set_ylim(0, 30)

        ax.plot(r/particle_diameter, G[t_index, :], marker='.')
        ax.text(0.7, 0.9, f'$\Delta t={t[t_index]}\mathrm{{s}}$', transform=ax.transAxes)
    

    common.save_gif(func, range(len(t)), fig, f'van_hove/figures_png/G_{file}.gif', fps=4)
    common.save_gif(func, range(len(t)), fig, f'/home/acarter/presentations/cin_first/figures/G_{file}.mp4', fps=4)


    