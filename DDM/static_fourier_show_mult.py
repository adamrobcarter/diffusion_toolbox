import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft
from DDM.static_fourier_show import show_static_fourier, form_factor_sphere

SAVE_TWO_D_PLOT = True
DISCRETE_COLORS = True

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    ymin = 1e100
    ymax = 0
    
    for i, file in enumerate(files := common.files_from_argv('DDM/data', 'static_fourier_')):
        if DISCRETE_COLORS:
            # color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:grey', 'tab:olive', 'tab:cyan'][i]
            color = f'C{i}'
        else:
            color = common.colormap(i, 0, len(files))

        particle_diameter, k_bin_mids, v = show_static_fourier(file, ax, color)

        ymin = np.nanmin([ymin, v[k_bin_mids>0.1].min()])
        ymax = np.nanmax([ymax, v[k_bin_mids>0.1].max()])
    
    ax.set_ylim(ymin*0.1, ymax*10)

    # common.add_exponential_index_indicator(ax, exponent=-4, anchor=(1, 1), xlabel='k')
    
    plt.legend(fontsize=9)

    filenames = '_'.join(files)
    filename = f'static_fourier_av_{filenames}'
    common.save_fig(fig, f'DDM/figures_png/{filename}.png', dpi=200)