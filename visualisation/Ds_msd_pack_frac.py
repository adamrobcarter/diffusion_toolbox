import matplotlib.pyplot as plt
import numpy as np
import common

fig, ax = plt.subplots(1, 1, figsize=(3.6, 3))

for file in common.files_from_argv('visualisation/data', 'Ds_from_MSD_first_'):
    try:
        data = common.load(f'visualisation/data/Ds_from_MSD_first_{file}.npz')
        Ds     = data['Ds']
        D_uncs = data['D_uncs']
        pack_frac = data['pack_frac_given']

        ax.errorbar(pack_frac, Ds, D_uncs, marker='o', linestyle='none', label=file)
    except Exception as err:
        print(err)
# ax.set_xticks(xs)

ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
ax.set_xlabel('packing fraction')
ax.legend(fontsize=7)

# common.save_fig(fig, '/home/acarter/presentations/cmd31/figures/Ds_hack.pdf', hide_metadata=True)
common.save_fig(fig, 'visualisation/figures_png/Ds_msd.png')

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