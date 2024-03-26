import matplotlib.pyplot as plt
import numpy as np
import common
import sys

collapse_x = True
collapse_y = True
collapse_x = False
collapse_y = False

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
    added_drift_x= data['added_drift_x']
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

    for box_size_index in range(N2_mean.shape[0]-1):
    # for L in [2**e for e in range(-2, 7)]:
        # L = box_sizes[box_size_index]

        delta_N_sq = N2_mean[box_size_index, :]
        t = np.copy(t_all)

        L = box_sizes[box_size_index]
        
        # computed theory interactions
        t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()))
        N2_theory = common.N2_nointer(t_theory, D0, N_var[box_size_index], L)
        
        if collapse_y:
            delta_N_sq /= N_var[box_size_index]
            N2_theory  /= N_var[box_size_index]
        if collapse_x:
            # t /= np.sqrt(Lx * Ly)
            t /= L**2
            t_theory /= L**2
            pass

        anomalous = delta_N_sq < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        if np.any(anomalous):
            print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
            delta_N_sq = delta_N_sq[~anomalous]
            t          = t         [~anomalous]

        label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
        # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

        
        exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=label, linestyle='none', marker='o', zorder=-1, markersize=5)
        # ax.plot(t_theory, N2_theory, color='black', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

    ax.legend(fontsize=7, loc='upper left')
    # ax.semilogy()
    # ax.semilogx()
    ax.set_xlabel('$t$')
    ax.set_ylabel('$\Delta N^2(t)$')
    title = f'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$'
    if not np.isnan(sigma):
        title += f', $\sigma={sigma:.3f}\mathrm{{\mu m}}$'
    if sigma_calced := data.get('particle_diameter_calced'):
        title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
    title += fr', $\nu_x={added_drift_x}\mathrm{{\mu m/s}}$'
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(f'box_counting/figures_png/N1N2_{file}.png', dpi=300)
