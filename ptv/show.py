import common
import numpy as np
import matplotlib.cm
import matplotlib.pyplot as plt

for file in common.files_from_argv('ptv/data/', 'ptv_'):
    data = common.load(f'ptv/data/ptv_{file}.npz')
    grid_xs = data['grid_xs']
    grid_ys = data['grid_ys']
    v_y       = data['v_y']
    v_x      = data['v_x']

    x_grid, y_grid = np.meshgrid(grid_xs, grid_ys, indexing='ij')

    fig, (ax_y, ax_x) = plt.subplots(2, 1, figsize=(5, 7.5))
    im = ax_y.pcolormesh(x_grid, y_grid, v_y, cmap=matplotlib.cm.seismic, shading='nearest',
                       vmin=-20, vmax=20)
    colorbar = fig.colorbar(im, ax=ax_y)
    colorbar.set_label('$v_y$, μm/s')

    im = ax_x.pcolormesh(x_grid, y_grid, v_x, cmap=matplotlib.cm.seismic, shading='nearest',
                       vmin=-20, vmax=20)
    colorbar = fig.colorbar(im, ax=ax_x)
    colorbar.set_label('$v_x$, μm/s')

    common.save_fig(fig, f'ptv/figures/ptv_y_{file}.png')
        