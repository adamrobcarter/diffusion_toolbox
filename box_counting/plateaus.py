import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_theory.nmsd, countoscope_theory.structure_factor

for file in common.files_from_argv('box_counting/data/', 'counted_'):
    fig, ax = plt.subplots(1, 1)
    fig2, ax2 = plt.subplots(1, 1)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean        = data['N2_mean']
    N2_std         = data['N2_std']
    phi            = data['pack_frac']
    sigma          = data['particle_diameter']
    time_step      = data['time_step']
    depth_of_field = data.get('depth_of_field')
    box_sizes      = data['box_sizes']
    N_mean         = data['N_mean']
    N_var          = data['N_var']
    N_mean_std     = data['N_mean_std']
    N_var_old      = data['N_var_old']
    N_var_std      = data['N_var_std']
    sep_sizes      = data['sep_sizes']


    num_timesteps = N2_mean.shape[1]
    plateaus = N2_mean[:, num_timesteps//2:-100].mean(axis=1)
    uncs     = N2_mean[:, num_timesteps//2:-100].std (axis=1)


    plot_plateau = ax.errorbar(sep_sizes[:20], plateaus[:20], yerr=0, zorder=10)
    plot_plateau = ax.errorbar(sep_sizes[:20], plateaus[:20], yerr=uncs[:20], color='tab:blue', alpha=0.7)
    # ax.set_xlim(sep_sizes.max()*1.1, sep_sizes.min()*1.1)
    ax.set_xlim(16, -16)
    ax.set_ylim(27, 35)
    common.save_fig(fig, f'/home/acarter/presentations/countoscope_may/plateaus_justhalf_{file}.png', dpi=200)

    plot_plateau = ax.errorbar(sep_sizes, plateaus, yerr=0, zorder=10, color='tab:blue')
    color_plateau = plot_plateau[0].get_color()
    plot_plateau = ax.errorbar(sep_sizes, plateaus, yerr=uncs, color='tab:blue', alpha=0.7)
    ax.tick_params(axis='y', colors=color_plateau)
    plat_ths = []
    for i in range(sep_sizes.size):
        plat_ths.append(
            countoscope_theory.nmsd.plateau_inter_2d(N_mean[i], box_sizes[i], lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)),
        )

    # ax.plot(sep_sizes, plat_ths, color='tab:red')

    ax.invert_xaxis()
    ax.set_xlabel('sep')
    ax.set_ylabel('plateau')

    ax.set_title(fr'$\phi={phi:.3f}$, $L={box_sizes.mean():.1f}\mathrm{{\mu m}}$')
    # ax.legend(loc='center left', fontsize=8)
    ax.set_xlim(16, -16)
    
    common.save_fig(fig, f'box_counting/figures_png/plateaus_just_{file}.png', dpi=200)
    common.save_fig(fig, f'/home/acarter/presentations/countoscope_may/plateaus_just_{file}.png', dpi=200)


    ax_mean = ax.twinx()
    ax_mean.plot([], [])
    plot_mean = ax_mean.errorbar(sep_sizes, N_mean, yerr=N_mean_std, alpha=0.3, color='tab:orange')
    plot_mean = ax_mean.errorbar(sep_sizes, N_mean, yerr=0, color='tab:orange')
    color_mean = plot_mean[0].get_color()
    ax_mean.tick_params(axis='y', colors=color_mean)

    ax_var = ax.twinx()
    ax_var.plot([], [], label='plateau')
    ax_var.plot([], [], label=r'$\langle N \rangle$')
    plot_var = ax_var.errorbar(sep_sizes, N_var, yerr=0,         color='tab:green', label='$\mathrm{Var}(N)$')
    plot_var = ax_var.errorbar(sep_sizes, N_var, yerr=N_var_std, color='tab:green', alpha=0.3)
    ax_var.plot(sep_sizes, N_var_old, color='tab:green', linestyle='dotted', label='Var old')
    color_var = plot_var[0].get_color()
    ax_var.tick_params(axis='y', colors=color_var)

    p = countoscope_theory.nmsd.plateau_inter_2d(N_mean.mean(), box_sizes.mean(), lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
    ax_var.hlines(p/2, sep_sizes.min(), sep_sizes.max(), color='tab:green', linestyle='dashed', label='var theory')

    ax_var.legend(loc='center left', fontsize=8)




    common.save_fig(fig, f'box_counting/figures_png/plateaus_{file}.png', dpi=200)
    common.save_fig(fig, f'/home/acarter/presentations/countoscope_may/plateaus_{file}.png', dpi=200)


    ax2.invert_xaxis()
    ax2.set_xlabel('sep')
    ax2.set_ylabel(r'plateau / $\langle N \rangle$')
    ax2.errorbar(sep_sizes, plateaus/N_mean, yerr=uncs/N_mean)
    common.save_fig(fig2, f'box_counting/figures_png/plateaus_rescaled_{file}.png', dpi=200)

    # print(data['countsa'].mean(axis=2))
    # print(data['countsa'].std(axis=2).mean())
    # print(N_mean[0])
    # print(data['countsb'].mean(axis=2))
    # print(data['countsb'].std(axis=2).mean())
    # print(N_mean[-1])