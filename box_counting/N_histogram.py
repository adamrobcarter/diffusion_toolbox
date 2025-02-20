import common
import matplotlib.pyplot as plt
import numpy as np
import visualisation.Ds_overlapped

LOSE_CORRELATED_SAMPLES = True

for file in common.files_from_argv('box_counting/data/', 'counted_counts_'):
    data = common.load(f'box_counting/data/counted_counts_{file}.npz')
    counts            = data['counts']
    box_sizes         = data['box_sizes']
    particle_diameter = data['particle_diameter']
    time_step         = data['time_step']

    num_rows = int(np.floor(np.sqrt(len(box_sizes))))
    num_cols = int(np.ceil(len(box_sizes) / num_rows))
    assert num_cols * num_cols >= len(box_sizes)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols*3, num_rows*3))
    axs_flat = axs.flatten()

    for box_size_index in range(len(box_sizes)):
        these_counts = counts[box_size_index, :, :, :]
        assert np.nanmax(these_counts) > 0

        if LOSE_CORRELATED_SAMPLES:
            L = box_sizes[box_size_index]
            D0, _, _ = visualisation.Ds_overlapped.get_D0(file)
            t_diffuseacrossbox = L**2 / (4 * D0)
            frames_diffuseacrossbox = int(np.ceil(t_diffuseacrossbox / time_step))
            these_counts = these_counts[:, :, ::frames_diffuseacrossbox]

        these_counts = these_counts.flatten()

        bins = np.arange(np.nanmin(these_counts)-0.5, np.nanmax(these_counts)+0.5+0.1, step=1) # +0.1 makes the arange inclusive of the end point not exclusive
        heights, _, _ = axs_flat[box_size_index].hist(these_counts, bins=bins)
        
        binmids = (bins[1:] + bins[:-1]) / 2
        hist_mean = np.average(binmids, weights=heights)
        hist_var  = np.average((binmids - hist_mean)**2, weights=heights)
        N_mean = np.nanmean(these_counts)
        N_var  = np.nanvar(these_counts)
        print(f'histogram mean = {hist_mean/N_mean} * N mean')
        print(f'histogram var  = {hist_var/N_var} * N var')

        if LOSE_CORRELATED_SAMPLES:
            axs_flat[box_size_index].set_title(f'$L={box_sizes[box_size_index]/particle_diameter:.3g}\sigma$ (nth t = {frames_diffuseacrossbox})')
        else:
            axs_flat[box_size_index].set_title(f'$L={box_sizes[box_size_index]/particle_diameter:.3g}\sigma$')
        axs_flat[box_size_index].set_xlabel('$N$')
        axs_flat[box_size_index].set_ylabel('count')

    fig.suptitle(file)
    filename = f'N_histogram_{file}'
    if LOSE_CORRELATED_SAMPLES:
        filename += '_losecorr'
    common.save_fig(fig, f'box_counting/figures_png/{filename}.png')