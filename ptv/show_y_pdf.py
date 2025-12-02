import common
import numpy as np
import matplotlib.cm
import matplotlib.pyplot as plt
import particle_detection.show
import preprocessing.frame1

def finite(arr):
    return arr[np.isfinite(arr)]

if __name__ == '__main__':
    for file in common.files_from_argv('ptv/data/', 'ptv_'):
        data = common.load(f'ptv/data/ptv_{file}.npz')
        grid_xs = data['grid_xs']
        grid_ys = data['grid_ys']
        v_x     = data['v_x']
        v_y     = data['v_y']
        n       = data['n'] # number of observations in each grid cell

        x_grid, y_grid = np.meshgrid(grid_xs, grid_ys, indexing='ij')
        
        # v_y = - np.nanmean(v_y, axis=0) # axis=0 means average over slices in the x?y? direction
        # v_y_unc = np.nanstd(v_y, axis=0) / np.sqrt(n.sum(axis=0))

        # # v_x = np.abs(v_x)
        # v_x = np.nanmean(v_x, axis=0)
        # v_x_unc = np.nanstd(v_x, axis=0) / np.sqrt(n.sum(axis=0))

        num_slices = v_y.shape[1]

        fig, (im_ax, y_ax, x_ax) = plt.subplots(1, 3, gridspec_kw={'width_ratios': (3, 1, 1)}, figsize=(11, 6), sharey='row')

        # particle_detection.show.go(file, im_ax, fig=None)
        preprocessing.frame1.go('faxtor006a_small', im_ax)

        v_y = - v_y
        assert np.nanmean(v_y) > 0, f'np.nanmean(v_y) = {np.nanmean(v_y)}'

        bins_y = np.linspace(-2, 15,   20)
        bins_x = np.linspace(-1.3,  1.3, 20)
        # bins_y = np.linspace(np.nanmin(v_y), np.nanmax(v_y), 20)
        # bins_x = np.linspace(np.nanmin(v_x), np.nanmax(v_x), 20)
        bin_y_centers = (bins_y[:-1] + bins_y[1:])/2
        bin_x_centers = (bins_x[:-1] + bins_x[1:])/2


        slice_width = data['window_size_x'] / num_slices

        for slice in range(num_slices):
            # print()
            print('slice', slice)

            data_y = finite(v_y[:, slice])
            data_x = finite(v_x[:, slice])

            if not data_y.size:
                print('no data y')
            if not data_x.size:
                print('no data x')

            hist_y, _ = np.histogram(data_y, bins=bins_y, density=True)
            hist_x, _ = np.histogram(data_x, bins=bins_x, density=True)
            # print('bins_y', bins_y)
            # print('v_y', finite(v_y[slice, :]))
            # print('hist', hist_y)

            # hist_y /= 5 # this is because with density=True, the area=1, not the max height
            hist_x /= 5 # but i want them all to be on the same scale, so can't just divide by the max

            offset = slice * slice_width
            y_ax.fill_between(bin_y_centers, offset, offset + hist_y*slice_width, color='tab:blue')
            x_ax.fill_between(bin_x_centers, offset, offset + hist_x*slice_width, color='tab:blue')

        y_ax.set_xlabel('$v_y$ (μm/s)')
        x_ax.set_xlabel('$v_x$ (μm/s)')

        im_ax.set_ylabel('$y$ (μm)')
        # ax1.set_ylabel('$y$ (μm)')

        y_ax.set_ylim(0, data['window_size_y'])
        x_ax.set_ylim(0, data['window_size_y'])

        y_ax.grid(axis='x')
        x_ax.grid(axis='x')

        y_ax.set_xlim(bin_y_centers.min(), bin_y_centers.max())
        x_ax.set_xlim(bin_x_centers.min(), bin_x_centers.max())
        # ax2.vlines(0, *ax2.get_ylim(), color='grey')

        print('v_SE', common.stokes_einstein_v(data['particle_diameter'], data['particle_material']))

        common.save_fig(fig, f'ptv/figures/ptv_y_pdf_{file}.png')
            