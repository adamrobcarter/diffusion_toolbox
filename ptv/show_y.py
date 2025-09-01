import common
import numpy as np
import matplotlib.cm
import matplotlib.pyplot as plt
import particle_detection.show
import preprocessing.frame1

for file in common.files_from_argv('ptv/data/', 'ptv_'):
    data = common.load(f'ptv/data/ptv_{file}.npz')
    grid_xs = data['grid_xs']
    grid_ys = data['grid_ys']
    v_x     = data['v_x']
    v_y     = data['v_y']
    n       = data['n'] # number of observations in each grid cell

    x_grid, y_grid = np.meshgrid(grid_xs, grid_ys, indexing='ij')
    
    v_y = - np.nanmean(v_y, axis=0) # axis=0 means average over slices in the x?y? direction
    v_y_unc = np.nanstd(v_y, axis=0) / np.sqrt(n.sum(axis=0))

    # v_x = np.abs(v_x)
    v_x = np.nanmean(v_x, axis=0)
    v_x_unc = np.nanstd(v_x, axis=0) / np.sqrt(n.sum(axis=0))

    fig, (im_ax, y_ax, x_ax) = plt.subplots(1, 3, gridspec_kw={'width_ratios': (3, 1, 1)}, figsize=(11, 6), sharey='row')

    # particle_detection.show.go(file, im_ax, fig=None)
    preprocessing.frame1.go('faxtor006a_small', im_ax)
    
    y_ax.errorbar(v_y, grid_ys, xerr=v_y_unc, linestyle='none', marker='o')
    x_ax.errorbar(v_x, grid_ys, xerr=v_x_unc, linestyle='none', marker='o')

    y_ax.set_xlabel('$v_y$ (μm/s)')
    x_ax.set_xlabel('$v_x$ (μm/s)')

    im_ax.set_ylabel('$y$ (μm)')
    # ax1.set_ylabel('$y$ (μm)')

    y_ax.set_ylim(0, data['window_size_y'])
    x_ax.set_ylim(0, data['window_size_y'])

    y_ax.grid(axis='x')
    x_ax.grid(axis='x')

    y_ax.set_xlim(0, 11)
    x_ax.set_xlim(-1.1, 1.1)
    # ax2.vlines(0, *ax2.get_ylim(), color='grey')

    print('v_SE', common.stokes_einstein_v(data['particle_diameter'], data['particle_material']))

    common.save_fig(fig, f'ptv/figures/ptv_y_{file}.png')
        