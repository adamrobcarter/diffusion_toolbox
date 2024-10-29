import common
import van_hove.van_hove as van_hove
import sys
import tqdm
import numpy as np
import matplotlib.pyplot as plt

for file in sys.argv[1:]:
    d = common.load(f"isf/data/F_{file}.npz")
    t     = d["t"]
    F     = d["F"]
    F_unc = d['F_unc']
    k     = d["k"]
    # particle_diameter = data.get('particle_diameter')
    
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    # ax.semilogy()

    def func(t_index):
        ax.clear()
        ax.set_ylabel('$F(k, t)$')
        ax.set_xlabel('$k$')
        # ax.semilogy() set_ylim(0.5, 30)
        ax.set_ylim(0, 1)
        ax.set_xlim(k[0, 0], k[0, -1])

        ax.plot(k[0, :], F[t_index, :]/F[0, :], marker='.', linestyle='none')
        ax.text(0.65, 0.9, f'$t={t[t_index]}\mathrm{{s}}$', transform=ax.transAxes)
    

    common.save_gif(func, range(len(t)), fig, f'isf/figures_png/F_{file}.gif', fps=4)
    # common.save_gif(func, range(len(t)), fig, f'/home/acarter/presentations/cin_first/figures/F_{file}.mp4', fps=4)


    