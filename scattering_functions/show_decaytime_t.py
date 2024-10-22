import common
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.optimize

NUM_Ks = 20


for file in common.files_from_argv('scattering_functions/data', 'F_'):


    fig, axs = plt.subplots(1, NUM_Ks, figsize=(3*NUM_Ks, 3.2))
    fig_single, ax_single = plt.subplots(1, 1)

    data = common.load(f'scattering_functions/data/F_{file}.npz')

    F_all = data['F']
    k_all = data['k']
    t_all = data['t']
    ks = k_all[0, :]
    sigma = data['particle_diameter']

    f_all = F_all / F_all[0, :]


    every_nth_k = int(math.ceil(k_all.shape[1] / NUM_Ks))
    every_nth_k = max(every_nth_k, 1)

    graph_i = 0

    displayed_ts = []
    
    for k_index in range(ks.size):
        display = False
        if k_index % every_nth_k == 0:
            display = True
        # if k_index < 10:
        #     display = True

        f = f_all[:, k_index]
        # t = t_all[k_index]
        k = ks[k_index]

        # if t > 1000:
        #     display = False

        decaytime = -t_all / np.log(f)

        if display:
            print(f'k={k:.2g}')
            ax = axs[graph_i]
            graph_i += 1

            good = decaytime > 0


            xs = t_all[good]
            ys = 1 / (k**2 * decaytime[good])
            ax.scatter(xs, ys)
            color = common.colormap(graph_i, 0, NUM_Ks)
            ax_single.plot(xs, ys, color=color)

            displayed_ts.append(k)

            if good.sum() == 0:
                continue

            x0 = xs[0]
            y0 = ys[0]

            # factor = 10
            # x1 = x0 * factor
            # y1_m1 = y0 * factor**(-1)
            # y1_m2 = y0 * factor**(-2)
            # y1_m3 = y0 * factor**(-4/3)

            # ax.plot([x0, x1], [y0, y1_m1], color='gray', linestyle='dotted')
            # ax.plot([x0, x1], [y0, y1_m2], color='gray', linestyle='dashed')
            # ax.plot([x0, x1], [y0, y1_m3], color='gray', linestyle='dashed')

            # early_start = 10
            # early_end = 20
            # late_start = 25
            # late_end = 35
            
            # func_early = lambda k, Lhydro_over_D0 : Lhydro_over_D0 / k
            # popt_early, pcov_early = scipy.optimize.curve_fit(func_early, k[good][early_start:early_end], decaytime[good][early_start:early_end])
            # k_th = k[good][early_start:early_end]
            # x_th = k_th * sigma
            # ax.plot(x_th, func_early(k_th, *popt_early), color=common.FIT_COLOR)

            # print(f'  e Lh/D0 = {popt_early[0]:.3f}')
            # Lh1 = 0.2 * sigma
            # Lh2 = 1 * sigma
            # print(f'  e D0 = {1/popt_early[0] * Lh1:.3f} (Lh={Lh1/sigma:.2g}σ)')
            # print(f'  e D0 = {1/popt_early[0] * Lh2:.3f} (Lh={Lh2/sigma:.2g}σ)')
            # target_D0 = 0.02
            # print(f'  for D0 = {target_D0}, need Lh={popt_early[0]*target_D0/sigma:.3g}σ')
            
            # if t < 200:
            #     func_late = lambda k, D0 : 1/ (D0 * k**2)
            #     popt_late, pcov_late = scipy.optimize.curve_fit(func_late, k[good][late_start:late_end], decaytime[good][late_start:late_end])
            #     k_th = k[good][late_start:late_end]
            #     x_th = k_th * sigma
            #     ax.plot(x_th, func_late(k_th, *popt_late), color=common.FIT_COLOR)

            #     print(f'  l D0 = {popt_late[0]:.3f}')


            ax.semilogx()
            ax.semilogy()
            ax.set_xlabel('$t$')
            ax.set_ylabel(r'$\tau=-t / \ln{F(k, t)}$')

            ax.set_title(f'$k={k:.1f}s$')
            # ax.set_ylim(0.5e0, 1e3)
        
    ax_single.semilogx()
    ax_single.semilogy()
    ax_single.set_xlabel('$t$')
    ax_single.set_ylabel(r'$\tau=-t / \ln{F(k, t)}$')
    # colorbar = common.colormap_colorbar(displayed_ts[0], displayed_ts[-1])
    colorbar = common.colormap_colorbar(0, NUM_Ks-1)
    fig_single.colorbar(colorbar, ax=ax_single, orientation='vertical', label='$k$ (s)?')



    common.save_fig(fig_single, f'scattering_functions/data/decaytime_single_t_{file}.png')
    common.save_fig(fig, f'scattering_functions/data/decaytime_t_{file}.png')