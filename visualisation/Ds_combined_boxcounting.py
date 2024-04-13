import common
import numpy as np
import matplotlib.pyplot as plt


fig, ax = plt.subplots(1, 1, figsize=(4, 4))

for file in (files := common.files_from_argv('visualisation/data', 'Ds_from_DDM_')):
    x = 0

    source_names = {
        'DDM': 'DDM',
        'f': '$f(k, t)$',
        'Fs': '$F_s(k, t)$',
        'boxcounting': 'Box Counting',
        'MSD': 'MSD'
    }

    all_Ds = []

    for source in ['boxcounting']:
        data = np.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
        Ds     = data['Ds']
        D_uncs = data['D_uncs']
        labels = data['labels']

        xs = x + np.linspace(-0.3, 0.3, num=Ds.size)

        ax.errorbar(xs, Ds, D_uncs, linestyle='none', marker='o', label=f'{file}')

        x += 1

        [all_Ds.append(D) for D in Ds]

    # ax.set_ylim(0, np.median(all_Ds)*2)
ax.set_ylabel('$D$')
ax.set_xticks([])
# fig.legend()
fig.legend(loc='upper right', bbox_to_anchor=(0.96, 0.9))
ax.set_title(f'errorbars not yet all correct')
filename = '_'.join(files)
common.save_fig(fig, f'visualisation/figures_png/Ds_combined_boxcounting_{filename}.png')