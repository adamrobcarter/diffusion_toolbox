import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_old as countoscope_theory.nmsd, countoscope_theory.structure_factor
import matplotlib.cm

for file in common.files_from_argv('box_counting/data/', 'counted_'):

    BOX_INDEX = -1

    fig, ax = plt.subplots(1, 1)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean        = data['N2_mean']
    N2_std         = data['N2_std']
    phi            = data['pack_frac']
    sigma          = data['particle_diameter']
    time_step      = data['time_step']
    depth_of_field = data.get('depth_of_field')

    sep_sizes = data['sep_sizes']

    all_counts = data['counts']

    counts = all_counts[BOX_INDEX, :, :,] # (x, y, t)
    counts = counts[:, ~np.isnan(counts).all(axis=(0, 2)), :] # remove all nan rows
    counts = counts[~np.isnan(counts).all(axis=(1, 2)), :, :] # remove all nan columns

    num_boxes_x = counts.shape[0]
    num_boxes_y = counts.shape[1]

    fig, axs = plt.subplots(num_boxes_x, num_boxes_y, figsize=(2*num_boxes_y, 2*num_boxes_x))

    max_y = 0

    for x in range(num_boxes_x):
        for y in range(num_boxes_y):
            # color = matplotlib.cm.afmhot(np.interp(counts[x, y, :].mean(), (counts.min(), counts.max()), (0.2, 0.75)))
            color = matplotlib.cm.afmhot(np.interp(counts[x, y, :].mean(), (np.percentile(counts, 25), np.percentile(counts, 75)), (0.2, 0.75)))
            axs[x][y].hist(counts[x, y, :], range=(counts.min(), counts.max()), color=color)
            max_y = max(max_y, axs[x][y].get_ylim()[1])

    # set all hists to have the same y-range     
    for x in range(num_boxes_x):
        for y in range(num_boxes_y):
            axs[x][y].set_ylim(0, max_y)
            axs[x][y].set_yticks([])
            axs[x][y].set_xticks([])

    common.save_fig(fig, f'box_counting/figures_png/N_hists_{file}.png', dpi=200)
