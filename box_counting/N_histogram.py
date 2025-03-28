import common
import matplotlib.pyplot as plt
import numpy as np
import visualisation.Ds_overlapped
import scipy.optimize
import scipy.stats

LOSE_CORRELATED_SAMPLES = False



def do_hist_fit(hist_x, hist_y, mean, var, ax=None):

    binmids = hist_x
    heights = hist_y
    binmids = binmids[np.isfinite(binmids)]
    heights = heights[np.isfinite(heights)]
    heights /= heights.sum() # normalise to 1

    hist_mean = np.average(binmids, weights=heights)
    hist_var  = np.average((binmids - hist_mean)**2, weights=heights)
    N_mean = mean
    N_var  = var 
        # print(f'histogram mean = {hist_mean/N_mean} * N mean')
        # print(f'histogram var  = {hist_var/N_var} * N var')

    if ax:
        ax.bar(binmids, heights, width=1)

        ax.set_xlabel('$N$')
        ax.set_ylabel('count')

    N_var_p0 = N_var * 1.5

    ##### unshifted
    fit_func = lambda k, l : scipy.stats.poisson.pmf(k, l)
    popt_unshifted, pcov_unshifted = scipy.optimize.curve_fit(fit_func, binmids, heights, p0=[N_mean])
    assert popt_unshifted[0] != N_mean
    if ax:
        ax.plot(binmids, fit_func(binmids, *popt_unshifted), label=f'Po(N, {popt_unshifted[0]:.3g})')

    ##### offset free parameter
    fit_func = lambda k, l, offset : scipy.stats.poisson.pmf(k, l, int(offset))
    popt_offsetfree, pcov_offsetfree = scipy.optimize.curve_fit(fit_func, binmids, heights, p0=(N_var*0.9, N_mean-N_var))
    
    print(f'  var_poisson/var fitshift', popt_offsetfree[0]/N_var)
    
    if ax:
        ax.plot(binmids, fit_func(binmids, *popt_offsetfree), label=f'Po(N-{popt_offsetfree[1]:.0f}, {popt_offsetfree[0]:.3g})', color='tab:pink')

    ##### offset locked
    fit_func = lambda k, l : scipy.stats.poisson.pmf(k, l, int(np.round(N_mean-N_var)))
    popt_offsetlocked, pcov_offsetlocked = scipy.optimize.curve_fit(fit_func, binmids, heights, p0=(N_var*0.9))
    print(f'  var_poisson/var nofitshift', popt_offsetlocked[0]/N_var)

    if ax:
        ax.plot(binmids, fit_func(binmids, *popt_offsetlocked), label=f'Po(N-{np.round(N_mean-N_var):.0f}, {popt_offsetlocked[0]:.3g})', color='tab:red')

        ax.legend(fontsize=8)

    return popt_offsetfree[0], 0

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'counted_'):
        data = common.load(f'box_counting/data/counted_{file}.npz')
        box_sizes         = data['box_sizes']
        particle_diameter = data['particle_diameter']
        time_step         = data['time_step']
        hist_x            = data['hist_x']
        hist_y            = data['hist_y']
        mean              = data['N_mean']
        var               = data['N_var']

        num_rows = int(np.floor(np.sqrt(len(box_sizes))))
        num_cols = int(np.ceil(len(box_sizes) / num_rows))
        assert num_cols * num_cols >= len(box_sizes)

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols*3, num_rows*3))
        axs_flat = axs.flatten()

        for box_size_index in range(len(box_sizes)):
            print(f'L = {box_sizes[box_size_index]/particle_diameter:.3g}sigma')

            axs_flat[box_size_index].set_title(f'$L={box_sizes[box_size_index]/particle_diameter:.3g}\sigma$')

            do_hist_fit(hist_x[box_size_index, :], hist_y[box_size_index, :],
                        mean[box_size_index], var[box_size_index], axs_flat[box_size_index])

        # fig.suptitle(file)
        filename = f'N_histogram_{file}'
        if LOSE_CORRELATED_SAMPLES:
            filename += '_losecorr'
        common.save_fig(fig, f'box_counting/figures_png/{filename}.png')