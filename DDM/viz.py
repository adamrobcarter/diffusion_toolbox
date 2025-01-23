import common
import matplotlib.pyplot as plt
import numpy as np

def show(file, d, t, k, time_origins):
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    k_indexes = np.array([43, 57, 71, 85]) / 100 * d.shape[2]
    for i, k_index in enumerate(k_indexes):
        k_index = int(k_index)
        ax = axs[i]
        extent = (t[1], t[-1], time_origins[-1], time_origins[0])
        to_plot = np.log(d[:, :, k_index])
        # d is shape (t0 * Delta t * k)
        vmin = np.nanquantile(to_plot, 0.05)
        vmax = np.nanquantile(to_plot, 0.95)
        # vmin = np.nanmax(to_plot)
        # vmax = np.nanmin(to_plot)
        # vmin = np.max(np.ma.masked_invalid(to_plot))
        # vmin = np.min(np.ma.masked_invalid(to_plot))
        # print(np.nanmax(to_plot), np.nanmin(to_plot), vmin, vmax)
        # im = ax.imshow(to_plot, vmin=vmin, vmax=vmax, interpolation='none', extent=extent)
        im = ax.pcolormesh(t, time_origins, to_plot)
        ax.set_ylim(ax.get_ylim()[1], ax.get_ylim()[0]) # flip y axis
        ax.set_xlabel('lag time $\Delta t$')
        ax.set_ylabel('real time $t$')
        ax.semilogx()
        ax.set_xlim(t[1], t[-1])
        # ax.set
        # ax.set_aspect(t[-1]/time_origins[-1])
        # ax.set_aspect(0.01)
        # ax.semilogx()
        um = f'\mathrm{{\mu m}}'
        ax.set_title(fr"$k={k[k_index]:.2f}{um}^{{-1}}$ ($\approx{2*np.pi/k[k_index]:.1f}{um}$)")
        fig.colorbar(im, ax=ax)
        fig.suptitle('$|I(q, t+\Delta t) - I(q, t)|^2$')

    common.save_fig(fig, f'DDM/figures_png/viz_{file}.png', dpi=300)

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data/', 'ddm_'):
        data = common.load(f'DDM/data/ddm_{file}.npz')
        d = data['F_D_sq_all'] # shape (real time) x (lag time) x (k)
        t = data['t']
        k = data['k']
        time_origins = data['time_origins']

        show(file, d, t, k, time_origins)
    