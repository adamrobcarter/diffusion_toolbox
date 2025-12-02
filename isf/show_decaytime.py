import common
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.optimize

DEFAULT_NUM_Ts = 10


def go(file, num_Ts, export_destination=None):

    fig, [axs] = plt.subplots(1, num_Ts, figsize=(3*num_Ts, 3.2), squeeze=False)
    fig_single, ax_single = plt.subplots(1, 1)

    data = common.load(f'isf/data/F_{file}.npz')

    F_all = data['F']
    k_all = data['k']
    t_all = data['t']
    k = k_all[0, :]
    sigma = data['particle_diameter']

    f_all = F_all / F_all[0, :]


    every_nth_k = int(math.ceil(k_all.shape[1] / num_Ts))
    every_nth_k = max(every_nth_k, 1)

    graph_i = 0

    displayed_ts = []
    
    for t_index in range(t_all.size):
        display = False
        if num_Ts > 1:
            if t_index % every_nth_k == 0:
                display = True
        else:
            if t_index == 1:
                display = True
        # if k_index < 10:
        #     display = True

        f = f_all[t_index, :]
        t = t_all[t_index]

        if t > 1000:
            display = False

        decaytime = -t / np.log(f)

        if display:
            print(f't={t}')
            ax = axs[graph_i]
            graph_i += 1

            good = decaytime > 0


            xs = k[good] * sigma
            ys = decaytime[good]
            ax.scatter(xs, ys)
            color = common.colormap(graph_i, 0, num_Ts)
            ax_single.plot(xs, ys, color=color)

            displayed_ts.append(t)

            if good.sum() == 0:
                continue

            # x0 = xs[0]
            # y0 = ys[0]

            # factor = 10
            # x1 = x0 * factor
            # y1_m1 = y0 * factor**(-1)
            # y1_m2 = y0 * factor**(-2)
            # y1_m3 = y0 * factor**(-4/3)

            # ax.plot([x0, x1], [y0, y1_m1], color='gray', linestyle='dotted')
            # ax.plot([x0, x1], [y0, y1_m2], color='gray', linestyle='dashed')
            # ax.plot([x0, x1], [y0, y1_m3], color='gray', linestyle='dashed')

            early_start = 11
            early_end = 22
            late_start = 26
            late_end = 35
            
            func_early = lambda k, Lhydro_over_D0 : Lhydro_over_D0 / k
            popt_early, pcov_early = scipy.optimize.curve_fit(func_early, k[good][early_start:early_end], decaytime[good][early_start:early_end])
            k_th = k[good][early_start:early_end]
            x_th = k_th * sigma
            ax.plot(x_th, func_early(k_th, *popt_early), color=common.FIT_COLOR)

            print(f'  early Lh/D0 = {popt_early[0]:.4f}')
            Lh1 = 0.2 * sigma
            Lh2 = 1 * sigma
            print(f'  early D0 = {1/popt_early[0] * Lh1:.4f} (Lh={Lh1/sigma:.2g}σ)')
            print(f'  early D0 = {1/popt_early[0] * Lh2:.4f} (Lh={Lh2/sigma:.2g}σ)')
            target_D0 = 0.02
            print(f'  early for D0 = {target_D0}, need Lh={popt_early[0]*target_D0/sigma:.3g}σ')
            
            if t < 200:
                func_late = lambda k, D0 : 1/ (D0 * k**2)
                popt_late, pcov_late = scipy.optimize.curve_fit(func_late, k[good][late_start:late_end], decaytime[good][late_start:late_end])
                k_th = k[good][late_start:late_end]
                x_th = k_th * sigma
                ax.plot(x_th, func_late(k_th, *popt_late), color=common.FIT_COLOR)

                print(f'  late  D0 = {popt_late[0]:.4f}')


            ax.semilogx()
            ax.semilogy()
            ax.set_xlabel(r'$k\sigma$')
            ax.set_ylabel(r'$\tau=-f / \ln{F(k, t)}$')

            ax.set_title(f'$t={t:.1f}s$')
            ax.set_ylim(0.2e0, 0.5e3)
        
    ax_single.semilogx()
    ax_single.semilogy()
    ax_single.set_xlabel(r'$k\sigma$')
    ax_single.set_ylabel(r'$\tau=-f / \ln{F(k, t)}$')
    # colorbar = common.colormap_colorbar(displayed_ts[0], displayed_ts[-1])
    colorbar = common.colormap_colorbar(0, num_Ts-1)
    fig_single.colorbar(colorbar, ax=ax_single, orientation='vertical', label='$t$ (s)')


    if export_destination:
        common.save_fig(fig, export_destination, hide_metadata=True)

    common.save_fig(fig_single, f'isf/data/decaytime_single_{file}.png')
    common.save_fig(fig, f'isf/data/decaytime_{file}.png')

if __name__ == '__main__':
    for file in common.files_from_argv('isf/data', 'F_'):
        go(file, DEFAULT_NUM_Ts)