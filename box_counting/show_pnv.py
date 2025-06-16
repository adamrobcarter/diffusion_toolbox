import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import particle_detection.add_drift_periodic

def go(file, N0_source='mean'):
    data = common.load(f'box_counting/data/pnv_{file}.npz')
    
    fig,       ax       = plt.subplots(1, 1)
    fig_slice, ax_slice = plt.subplots(1, 1)

    used_avs = []
    used_Ls = []

    for box_size_index in range(len(data['box_sizes_x'])):
        # box_size_index = 0
        N1N2 = data['N1N2_per_box'][box_size_index, :, :, :]
        x_grid = data['box_coords_1'][box_size_index, :, :, 0]
        y_grid = data['box_coords_1'][box_size_index, :, :, 1]
        y = y_grid[0, :]
        x = x_grid[:, 0]

        # print('N1N2 avg', np.nanmean(N1N2))

        # some of the above can be nan (because different box_size_indexes can have different numbers of boxes)
        # so lets remove nan values
        finite_x = np.isfinite(x)
        finite_y = np.isfinite(y)
        x = x[finite_x]
        y = y[finite_y]
        # print(N1N2.shape)
        # print(finite_x.shape, finite_y.shape)
        # print(np.ix_(finite_x, finite_y))
        # print(N1N2[np.ix_(finite_x, finite_y)].shape)
        N1N2 = N1N2[np.ix_(finite_x, finite_y)] # this is really N1N2[np.ix_(finite_x, finite_y), :]
        # N1N2 = N1N2[finite_x, finite_y, :]
        # print(N1N2.shape)
        x_grid = x_grid[np.ix_(finite_x, finite_y)]
        y_grid = y_grid[np.ix_(finite_x, finite_y)]

        L1 = data['box_sizes_x'][box_size_index]
        L2 = data['box_sizes_y'][box_size_index]
        dt = data['time_step']
        N_mean = data['N_mean'][box_size_index]
        N_var = data['N_var'][box_size_index]

        # N0 = N_mean * (N_mean/N_var)
        # N0 = N_var * 1.8
        if N0_source == 'var':
            N0 = N_var
        elif N0_source == 'mean':
            N0 = N_mean
        elif N0_source == 'special':
            N0 = N_var**2 / N_mean
        elif N0_source == 'N0S0':
            N0 = N_mean * common.S_k_zero(data['pack_frac'])
        # N0 /= common.S_k_zero(data['pack_frac'])
        assert np.all(N1N2[:, : 0] == 0)

        # print('mean/var', data['N_mean'][box_size_index]/data['N_var'][box_size_index])
        v_o1 = L1 * N1N2[:, :, 1] / (N0 * dt)
        v_o2 = L1 * N1N2[:, :, 1] / (N0 * dt) / (1 - np.sqrt(4*0.04*dt)/L2/np.sqrt(np.pi))
        v_o3 = L1 * N1N2[:, :, 1] / (N0 * dt) / (1 - np.sqrt(4*0.04*dt)/L2/np.sqrt(np.pi) - 0)
        
        # what is this 70 thing? wait it's just for that visualisation right?
        slice = np.index_exp[:, 70]
        if 'v_profile' in data:
            slice = np.index_exp[70, :]
        
        v_o1_slice = np.copy(v_o1[slice])
        v_o2_slice = np.copy(v_o2[slice])
        v_o1_slice = v_o1.mean(axis=0) # is this mean over time?
        v_o2_slice = v_o2.mean(axis=0)
        v_o1[slice] = -np.abs(v_o1).max()

        if True:
            used_avs.append(v_o2.mean())
            used_Ls.append(L1)

        min_L_theory = 0

        if 'v_profile' in data:
            # print(y)
            # v_on_slice = np.cos(2*np.pi * y / data['window_size_y'])
            v_func = particle_detection.add_drift_periodic.get_v[str(data['v_profile'])]
            v_theory_on_slice = v_func(data['window_size_x'], data['window_size_y'], x_grid, y_grid)[slice]
            # print(v_on_slice)

            ax_slice.plot(y, v_theory_on_slice[:, 0], color='black', linestyle='dashed')

            # ratio  = v_on_slice[:, 0] / v_o1_slice
            # ratio2 = v_on_slice[:, 0] / v_o2_slice
            # print('ratios', np.nanmean(ratio[ratio > 0]), np.nanmean(ratio2[ratio > 0]))
            # print('1/ratios', 1/np.nanmean(ratio[ratio > 0]), 1/np.nanmean(ratio2[ratio > 0]))

            v_theory_max = np.abs(v_theory_on_slice).max()
            min_L_theory = v_theory_max * dt
            # print('min L theory', min_L_theory)

        # print(x_grid.shape, y_grid.shape, v_o1.shape)
        ax.pcolormesh(x_grid, y_grid, v_o1, cmap=matplotlib.cm.seismic, vmin=-np.abs(v_o1).max(), vmax=np.abs(v_o1).max(), shading='nearest')

        # plot real data
        alpha = 1.0 if L1 > min_L_theory else 0.3
        # print(L1 > min_L_theory, alpha)
        ax_slice.plot(y, v_o1_slice, label=f'L={L1:.1f}', color=common.colormap(box_size_index, 0, len(data['box_sizes_x'])),
                      alpha=alpha)#'order 1')
        # ax_slice.plot(y, v_o2_slice, label='order 2')

    ax_slice.hlines(0, *ax.get_xlim(), color='black', linestyle=':')
    ax_slice.set_ylabel(r'$v$ ($\mathrm{\mu m / s}$)')

    ax_slice.legend()
    common.save_fig(fig, f'box_counting/figures_png/pnv_{file}.png')
    common.save_fig(fig_slice, f'box_counting/figures_png/pnv_slice_{file}.png')

    print(f'<N1N2 - N2N1> = {common.format_val_and_unc(np.mean(used_avs), np.std(used_avs))}')
    return used_avs, used_Ls, dt

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'pnv_'):
        go(file)