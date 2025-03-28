import common
import matplotlib.pyplot as plt
import numpy as np
import visualisation.Ds_overlapped
import scipy.optimize
import scipy.stats

LOSE_CORRELATED_SAMPLES = False

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
        print()
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
        heights, _, _ = axs_flat[box_size_index].hist(these_counts, bins=bins, density=True)
        
        binmids = (bins[1:] + bins[:-1]) / 2
        hist_mean = np.average(binmids, weights=heights)
        hist_var  = np.average((binmids - hist_mean)**2, weights=heights)
        N_mean = np.nanmean(these_counts)
        N_var  = np.nanvar(these_counts)
        # print(f'histogram mean = {hist_mean/N_mean} * N mean')
        # print(f'histogram var  = {hist_var/N_var} * N var')

        if LOSE_CORRELATED_SAMPLES:
            axs_flat[box_size_index].set_title(f'$L={box_sizes[box_size_index]/particle_diameter:.3g}\sigma$ (nth t = {frames_diffuseacrossbox})')
        else:
            axs_flat[box_size_index].set_title(f'$L={box_sizes[box_size_index]/particle_diameter:.3g}\sigma$')
        axs_flat[box_size_index].set_xlabel('$N$')
        axs_flat[box_size_index].set_ylabel('count')

        fit_func = lambda k, l : scipy.stats.poisson.pmf(k, l)
        popt, pcov = scipy.optimize.curve_fit(fit_func, binmids, heights, p0=[N_mean*0.9])
        axs_flat[box_size_index].plot(binmids, fit_func(binmids, *popt), label=f'Po(N, {popt[0]:.3g})')
        print()

        fit_func = lambda k, l, offset : scipy.stats.poisson.pmf(k, l, int(offset))
        popt, pcov = scipy.optimize.curve_fit(fit_func, binmids, heights, p0=(N_var*0.9, N_mean-N_var))
        # # print(f'var = {common.format_val_and_unc(popt[0], np.sqrt(pcov[0,0]))}')
        print(f'var_poisson/var fitshift', popt[0]/N_var)
        print(pcov)
        axs_flat[box_size_index].plot(binmids, fit_func(binmids, *popt), label=f'Po(N-{popt[1]:.0f}, {popt[0]:.3g})', color='tab:red')
        # print('N_mean-N_var', N_mean-N_var)
        # print('N_mean N_var', N_mean, N_var)
        # print('peak', fit_func(N_mean, N_var, N_mean-N_var))
        # axs_flat[box_size_index].plot(binmids, fit_func(binmids, N_mean, 0), label='Nvar3')
        # axs_flat[box_size_index].plot(binmids, fit_func(binmids, N_var, N_mean-N_var), label='Nvar2')
        # # print(binmids, fit_func(binmids, N_var))
        # axs_flat[box_size_index].legend()
        print()

        

        fit_func = lambda k, l : scipy.stats.poisson.pmf(k, l, int(np.round(N_mean-N_var)))
        popt, pcov = scipy.optimize.curve_fit(fit_func, binmids, heights, p0=(N_var*0.9))
        print(f'var_poisson/var nofitshift', popt[0]/N_var)
        print(pcov)
        axs_flat[box_size_index].plot(binmids, fit_func(binmids, *popt), label=f'Po(N-{np.round(N_mean-N_var):.0f}, {popt[0]:.3g})', color='tab:red')

        axs_flat[box_size_index].legend(fontsize=8)
        print()

    # fig.suptitle(file)
    filename = f'N_histogram_{file}'
    if LOSE_CORRELATED_SAMPLES:
        filename += '_losecorr'
    common.save_fig(fig, f'box_counting/figures_png/{filename}.png')

    print('TODO: YOU NO LONGER NEED TO USE COUNTED_COUNTS now THE HISTS ARE IN NORMAL')