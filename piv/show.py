import common
import numpy as np
import matplotlib.cm
import matplotlib.pyplot as plt
import tqdm

for file in common.files_from_argv('piv/data/', 'ptv_'):
    data = common.load(f'piv/data/piv_{file}.npz')
    # grid_xs = data['grid_xs']
    # grid_ys = data['grid_ys']
    v_y      = -data['v_y']
    v_x      =  data['v_x']
    signal_to_noise = data['signal_to_noise']

    # import openpiv.validation
    # for frame in tqdm.trange(v_y.shape[0], leave=False):
    #     invalid = openpiv.validation.sig2noise_val(signal_to_noise[frame], threshold = 1.05)

    #     # if frame == 0:
    #     #     print('sig to noise')
    #     #     common.term_hist(signal_to_noise[frame].flatten(), bins=np.linspace(0, 100))
    #     #     print('invalid', invalid.sum()/invalid.size)

    #     assert invalid.sum() < invalid.size
    #     v_x[frame, invalid] = np.nan
    #     v_y[frame, invalid] = np.nan
    
    assert np.isfinite(v_x).any()
    assert np.isfinite(v_y).any()
    v_x = np.nanmean(v_x, axis=0)
    v_y = np.nanmean(v_y, axis=0)
    # this averaging is wrong because if one frame gives a velocity and the next has no particle so
    # gives zero, it will average to v/2

    # x_grid, y_grid = np.meshgrid(grid_xs, grid_ys, indexing='ij')

    fig, (ax_y, ax_x) = plt.subplots(2, 1, figsize=(5, 7.5))
    # im = ax_y.pcolormesh(x_grid, y_grid, v_y, cmap=matplotlib.cm.seismic, shading='nearest',
    #                    vmin=-20, vmax=20)
    im = ax_y.imshow(v_y, cmap=matplotlib.cm.seismic, vmin=-20, vmax=20)
    colorbar = fig.colorbar(im, ax=ax_y)
    colorbar.set_label('$v_y$, μm/s')

    # im = ax_x.pcolormesh(x_grid, y_grid, v_x, cmap=matplotlib.cm.seismic, shading='nearest',
    #                    vmin=-20, vmax=20)
    im = ax_x.imshow(v_x, cmap=matplotlib.cm.seismic, vmin=-20, vmax=20)
    colorbar = fig.colorbar(im, ax=ax_x)
    colorbar.set_label('$v_x$, μm/s')

    common.save_fig(fig, f'piv/figures/piv_{file}.png')
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    v_mag = np.sqrt(v_x**2 + v_y**2)
    mask = v_mag < 5
    v_y[mask] = np.nan
    v_x[mask] = np.nan
    ax.quiver(v_x, v_y, np.arctan2(v_x, v_y), angles='xy', scale_units='xy', scale=10)
    ax.set_aspect('equal')
        
    common.save_fig(fig, f'piv/figures/piv_arrows_{file}.png', dpi=600)