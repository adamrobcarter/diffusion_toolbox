import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt

END_TIME_FRAC = 0.6

def go(file):
    data = common.load(f'van_hove/data/G_{file}.npz')
    G = data['G']
    r = data['r']
    t = data['t']
    particle_diameter = data.get('particle_diameter')
    
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    # ax.semilogy()

    def func(t_index):
        ax.clear()
        ax.set_ylabel('$G(r, t)$')
        ax.set_xlabel('$r/\sigma$')
        # ax.semilogy() set_ylim(0.5, 30)
        ax.set_ylim(0.4, 70)
        ax.semilogy()

        ax.plot(r/particle_diameter, G[t_index, :], marker='.', color='tab:green')
        ax.text(0.65, 0.9, f'$t={t[t_index]}\mathrm{{s}}$', transform=ax.transAxes)
    

    common.save_gif(func, range(int(len(t)*END_TIME_FRAC)), fig, f'van_hove/figures_png/G_{file}.gif', fps=4)
    common.save_gif(func, range(int(len(t)*END_TIME_FRAC)), fig, f'/home/acarter/presentations/cmd31/figures/G_{file}.mp4', fps=4)


    

if __name__ == '__main__':
    for file in sys.argv[1:]:
        go(file)