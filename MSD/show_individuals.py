import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize
import matplotlib.cm

for file in common.files_from_argv('MSD/data', 'msd_individuals_'):
    data = np.load(f'MSD/data/msd_individuals_{file}.npz')
    msds = data['msds']

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3))

    t = np.arange(0, msds.shape[1]) * data['time_step']

    # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
    for i in range(100):
        ax.plot(t[1:], msds[i, 1:], marker='.', markersize=3, linestyle='none')


    # ax.loglog()
    # ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$r(t)^2$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$\Delta t$ (s)')

    common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/msd_individuals_{file}.png', dpi=200, hide_metadata=True)
    # we save as png not pdf cause with all those lines the pdf will take a while to render
    common.save_fig(fig, f'MSD/figures_png/msd_individuals_{file}.png')