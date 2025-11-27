import matplotlib.pyplot as plt
import numpy as np
import common
import matplotlib.cm
import scipy.optimize
import countoscope_theory.n1n2
import countoscope_theory.structure_factor
from visualisation.Ds_overlapped import get_D0

# collapse_x = False
# collapse_y = False

COLLAPSE_Y_NONE = 0
COLLAPSE_Y_N0S0 = 1
COLLAPSE_Y_N0 = 2

collapse_x = True
collapse_y = COLLAPSE_Y_N0
# collapse_y = COLLAPSE_Y_NONE

SHOW_SHORT_TIME_EXPANSION = False
SHOW_FULL_THEORY = True

ZOOM = True

def go(file):
    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/pnv_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N1N2_mean = data['N1N2']
    N1N2_std  = data['N1N2_std']
    phi       = data['pack_frac']
    sigma     = data['particle_diameter']
    # drift_x   = data.get('drift_x', 0)
    drift_x = data['velocity_multiplier']
    # drift_y   = data.get('drift_y', 0)
    box_sizes = data['box_sizes_x']
    time_step = data['time_step']
    N_mean    = data['N_mean']
    N_var     = data['N_var']


    num_timesteps = N1N2_mean.shape[1]
    num_boxes     = N1N2_mean.shape[0]
    t_all = np.arange(0, num_timesteps) * time_step

    # reduce = 1
    # t        = t      [::reduce]
    # # t_theory = np.logspace()
    # N2_mean  = N2_mean[:, ::reduce]
    # N2_std   = N2_std [:, ::reduce]

    # N2_mean = N2_mean[:, :N2_mean.shape[1]//2]
    # t_all = t_all[:N2_mean.shape[1]]

    D_MSD, _, _ = get_D0(file.split('_drifted')[0])

    for box_size_index in range(N1N2_mean.shape[0]):
    # for L in [2**e for e in range(-2, 7)]:
        # L = box_sizes[box_size_index]

        # if box_size_index % 3 != 0:
        #     continue

        N1N2 = N1N2_mean[box_size_index, :]
        t = np.copy(t_all)

        L = box_sizes[box_size_index]
        
        # computed theory interactions
        # t_theory = np.copy(t)[:10:3]
        t_vx_over_L_theory = np.linspace(0, 0.25, num=10)
        t_theory = t_vx_over_L_theory * L / drift_x
        
        N_theory_shorttime = N_mean[box_size_index]*drift_x*t_theory/L # Grace's presentation slide 7

        # fit
        # func = lambda t, drift, D: common.N1N2_square(t, D, N_mean[box_size_index], L, drift)
        # popt, pcov = scipy.optimize.curve_fit(func, t, N1N2)
        # print(popt)

        if SHOW_FULL_THEORY:
            # full_theory = countoscope_theory.n1n2.N1N2_inter_old(t_theory, D_MSD, N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma), drift_x)
            full_theory = countoscope_theory.n1n2.N1N2_inter(t_theory, D_MSD, N_mean[box_size_index], L, phi, sigma, drift_x)
            assert np.isfinite(full_theory).all(), f'Non-finite values in theory: {full_theory}'

        xlabel = '$t$'
        ylabel = '$N_2(t)N_1(0) - N_1(0)N_2(t)$'
    
        if collapse_y == COLLAPSE_Y_N0:
            rescale_y = N_mean[box_size_index]
            ylabel = '$[N_2(t)N_1(0) - N_1(0)N_2(t)] / N_0$'

        elif collapse_y == COLLAPSE_Y_N0S0:
            # rescale_y = L**2
            rescale_y = N_mean[box_size_index] * common.S_k_zero(phi)
            ylabel = '$[N_2(t)N_1(0) - N_1(0)N_2(t)] / N_0S_0$'

        else:
            rescale_y = 1
            ylabel = '$[N_2(t)N_1(0) - N_1(0)N_2(t)]$'

        N1N2 /= rescale_y
        # N_theory   /= rescale_y
        N_theory_shorttime /= rescale_y
        if SHOW_FULL_THEORY:
            full_theory /= rescale_y
    
        if collapse_x and drift_x > 0:
            rescale_x = drift_x / L

            # t /= np.sqrt(Lx * Ly)
            t *= rescale_x
            t_theory *= rescale_x
            xlabel = r'$t \nu_x/L$'

            # ax.set_xlim(t[0], t[-1])

        # label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
        label = rf'$L={L/sigma:.1f}\sigma$'
        # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

        color = common.colormap(box_size_index, 0, len(box_sizes))  
        ax.plot(t, N1N2, label=label, linestyle='-', marker='o', zorder=-1, markersize=5, color=color)
        
        n = np.arange(1, N1N2_std.shape[1]+1)
        n = n[::-1]
        err = N1N2_std[box_size_index, :] / np.sqrt(n)
        # ax.fill_between(t, N1N2-err, N1N2+err, color=color, alpha=0.2)

        # ax.plot(t_theory, N_theory, color=common.FIT_COLOR, linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        if SHOW_SHORT_TIME_EXPANSION:
            ax.plot(t_theory, N_theory_shorttime, linestyle='dotted', color='black', label=r'short time exp. $\langle N \rangle$' if box_size_index==0 else None)

        if SHOW_FULL_THEORY:
            ax.plot(t_theory, full_theory, linestyle='dashed', color='black', label=r'theory (interac.) num. int.' if box_size_index==0 else None)

    ax.legend(fontsize=7)
    # ax.semilogy()
    # ax.semilogx()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(fr'$\phi = {phi:.1g}$')

    if drift_x != 0:
        # ax.set_xlim(0, 4)
        pass

    # if 
    # ax.set_ylim(-0.0005, 0.0025)
    if collapse_x:
        ax.set_xlim(-0.02, 2)

    if collapse_y == COLLAPSE_Y_N0S0:
        ax.set_ylim(-0.01, 0.5)
    elif collapse_y == COLLAPSE_Y_N0:
        ax.set_ylim(-0.01, 0.6)
    else:
        ax.set_ylim(-0.01, 0.2)

    if ZOOM:
        ax.set_xlim(-0.01, 0.3)
        ax.set_ylim(-0.02, 0.2)

    title = fr'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$'
    if not np.isnan(sigma):
        title += fr', $\sigma={sigma:.3f}\mathrm{{\mu m}}$'
    if sigma_calced := data.get('particle_diameter_calced'):
        title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
    title += fr', $\nu_x={drift_x}\mathrm{{\mu m/s}}$'
    # ax.set_title(title)

    fig.tight_layout()
    common.save_fig(fig, f'box_counting/figures_png/N1N2_{file}.png', dpi=300)


if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data', 'pnv_'):
        go(file)