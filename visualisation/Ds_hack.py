import matplotlib.pyplot as plt
import numpy as np
import common

fig, ax = plt.subplots(1, 1, figsize=(8, 3))

xs = np.array([0.02, 0.34, 0.66])
msd = np.array([4.34, 3.16, 1.85])
msd_uncs = np.array([0.13, 0.09, 0.06])
n = np.array([4.16, 3.10, 1.75])
n_uncs = np.array([0.28, 0.22, 0.11])

ax.errorbar(xs+0.005, n, n_uncs, marker='o', linestyle='none', label='Countoscope', color='tab:orange')
ax.errorbar(xs-0.005, msd, msd_uncs, marker='o', linestyle='none', label='MSD', color='tab:blue')
ax.set_xticks(xs)

ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
ax.set_xlabel('packing fraction')
ax.legend()

common.save_fig(fig, '/home/acarter/presentations/cmd31/figures/Ds_hack.pdf', hide_metadata=True)
common.save_fig(fig, 'visualisation/figures_png/Ds_hack.png')

"""

xs = np.array([0.02, 0.34, 0.66])
msd = np.array([4.34, 3.16, 1.85])
msd_uncs = np.array([0.13, 0.09, 0.06])
n = np.array([4.16, 3.10, 1.75])
n_uncs = np.array([0.28, 0.22, 0.11])

for p in [0, 1, 2]:
    fig, ax = plt.subplots(1, 1, figsize=(3.2, 3))


    ax.errorbar(xs[p]-0.02, msd[p], msd_uncs[p], marker='o', linestyle='none', label='MSD')
    ax.errorbar(xs[p]+0.02, n[p], n_uncs[p], marker='o', linestyle='none', label='Countoscope')
    # ax.set_xticks(xs)

    ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
    # ax.set_xlabel('packing fraction')
    ax.legend()

    common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/Ds_hack_{p}.pdf', hide_metadata=True)
    common.save_fig(fig, f'visualisation/figures_png/Ds_hack_{p}.png')
    """