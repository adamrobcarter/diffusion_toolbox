import common
import numpy as np
import matplotlib.cm
import matplotlib.pyplot as plt

for file in common.files_from_argv('ptv/data/', 'ptv_'):
    data = common.load(f'ptv/data/ptv_{file}.npz')
    grid_xs = data['grid_xs']
    grid_ys = data['grid_ys']
    v       = data['v']

    x_grid, y_grid = np.meshgrid(grid_xs, grid_ys, indexing='ij')

    fig, ax = plt.subplots(1, 1)
    im = ax.pcolormesh(x_grid, y_grid, v, cmap=matplotlib.cm.seismic, shading='nearest',
                       vmin=-20, vmax=20)
    colorbar = fig.colorbar(im)
    colorbar.set_label('$v$, px/frame')

    common.save_fig(fig, f'ptv/figures/ptv_{file}.png')
        