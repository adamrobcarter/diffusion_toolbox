import matplotlib.pyplot as plt
import numpy as np
import common
import sys

collapse_x = False
collapse_y = False
collapse_x = True
collapse_y = True

for file in sys.argv[1:]:

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/countedN1N2_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean      = data['N2_mean']
    N2_std       = data['N2_std']
    N_stats      = data['N_stats']
    phi          = data['pack_frac']
    sigma        = data['particle_diameter']
    drift_x      = data.get('drift_x', 0)
    drift_y      = data.get('drift_y', 0)
    box_sizes = data['box_sizes']
    time_step    = data['time_step']
    N_mean    = N_stats[:, 1]
    N_var     = N_stats[:, 2]

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t_all = np.arange(0, num_timesteps) * time_step

    # reduce = 1
    # t        = t      [::reduce]
    # # t_theory = np.logspace()
    # N2_mean  = N2_mean[:, ::reduce]
    # N2_std   = N2_std [:, ::reduce]

    # N2_mean = N2_mean[:, :N2_mean.shape[1]//2]
    # t_all = t_all[:N2_mean.shape[1]]

    # D0 = { # countoscope paper, table 1
    #     'alice0.02': 0.0416,
    #     'alice0.02_overlapped': 0.0416,
    #     'alice0.34': 0.0310,
    #     'alice0.66': 0.0175
    # }[file]
    D0 = 0.0416

    for box_size_index in range(N2_mean.shape[0]):
    # for L in [2**e for e in range(-2, 7)]:
        # L = box_sizes[box_size_index]

        N1N2 = N2_mean[box_size_index, :]
        t = np.copy(t_all)

        L = box_sizes[box_size_index]
        
        # computed theory interactions
        t_theory = np.copy(t)
        # N_theory = common.N1N2_nointer(t_theory, D0, N_var[box_size_index], (L, L), (added_drift_x, added_drift_y))
        N_theory = common.N1N2_square(t_theory, D0, N_mean[box_size_index], L, drift_x)
        # print(t_theory)
        # print(N_theory)
        N_theory_shorttime = N_mean[box_size_index]*drift_x*t_theory/L # Grace's presentation slide 7

        if collapse_y:
            N1N2 /= L**2
            N_theory   /= L**2
            N_theory_shorttime /= L**2
        if collapse_x and drift_x > 0:
            # t /= np.sqrt(Lx * Ly)
            t *= drift_x / L
            t_theory *= drift_x / L

            ax.set_xlim(t[0], t[-1])

        anomalous = N1N2 < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        if np.any(anomalous):
            print(f'found {anomalous.sum()/N1N2.size*100:.3f}% anomalous')
            # delta_N_sq = delta_N_sq[~anomalous]
            # t          = t         [~anomalous]

        # label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
        label = rf'$L={L/sigma:.1f}\sigma$'
        # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

        
        exp_plot = ax.plot(t, N1N2, label=label, linestyle='none', marker='o', zorder=-1, markersize=5)
        ax.plot(t_theory, N_theory, color='black', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        # ax.plot(t_theory, N_theory_shorttime, linestyle='dotted', color='black')

    ax.legend(fontsize=7, loc='upper left')
    # ax.semilogy()
    # ax.semilogx()
    xlabel = r'$t \nu_x/L$' if collapse_x else '$t$'
    ylabel = '$[N_1(t)N_2(0) - N_1(0)N_2(t)] / L^2$' if collapse_y else '$N_1(t)N_2(0) - N_1(0)N_2(t)$'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    title = f'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$'
    if not np.isnan(sigma):
        title += f', $\sigma={sigma:.3f}\mathrm{{\mu m}}$'
    if sigma_calced := data.get('particle_diameter_calced'):
        title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
    title += fr', $\nu_x={drift_x}\mathrm{{\mu m/s}}$'
    ax.set_title(title)

    fig.tight_layout()
    common.save_fig(fig, f'box_counting/figures_png/N1N2_{file}.png', dpi=300)
