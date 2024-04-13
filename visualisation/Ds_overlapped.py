import common
import numpy as np
import matplotlib.pyplot as plt

for file in common.files_from_argv('visualisation/data', 'Ds_from_DDM_'):
    x = 0
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    source_names = {
        'DDM': 'DDM',
        'f': '$f(k, t)$',
        'Fs': '$F_s(k, t)$',
        'boxcounting': 'Box Counting',
        'MSD': 'MSD'
    }

    all_Ds = []

    for source in ['DDM', 'f', 'Fs', 'boxcounting', 'MSD']:
        data = np.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
        Ds     = data['Ds']
        D_uncs = data['D_uncs']

        if source in ['f', 'Fs', 'DDM']:
            xs = data['ks']
        if source == 'boxcounting':
            xs = 2 * np.pi / data['Ls']
        if source == 'MSD':
            xs = 1

        ax.errorbar(xs, Ds, D_uncs, linestyle='none', marker='o', label=f'$D$ from {source_names[source]}')

        [all_Ds.append(D) for D in Ds]

    # ax.set_ylim(0, np.median(all_Ds)*2)
    ax.set_ylabel('$D$')
    ax.set_xticks([])
    ax.semilogx()
    # fig.legend()
    fig.legend(loc='upper right', bbox_to_anchor=(0.96, 0.9))
    ax.set_title(f'{file}, errorbars not yet all correct')
    common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_{file}.png')