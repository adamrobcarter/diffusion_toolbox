import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
import sys

print('only use me for low densities!')

for file in sys.argv[1:]:

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean      = data['N2_mean']
    N2_std       = data['N2_std']
    N_stats      = data['N_stats']
    phi          = data['pack_frac']
    sigma        = data['particle_diameter']
    sigma_calced = data['particle_diameter_calced']
    time_step    = data['time_step']

    box_sizes = N_stats[:, 0]
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

    t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()))

    # D0 = { # countoscope paper, table 1
    #     0.02: 0.0416,
    #     0.34: 0.0310,
    #     0.66: 0.0175
    # }[phi]

    for box_size_index, L in enumerate(box_sizes):
    # for L in [2**e for e in range(-2, 7)]:
        L = box_sizes[box_size_index]

        delta_N_sq = N2_mean[box_size_index, :]
        t = np.copy(t_all)

        anomalous = delta_N_sq < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        if np.any(anomalous):
            print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
            delta_N_sq = delta_N_sq[~anomalous]
            t          = t         [~anomalous]


        # t = np.arange(0, len(delta_N_sq))[1:]/2
        # delta_N_sq = delta_N_sq # [1:] is because the point at t=0 msd=0 plots weirdly

        
        L_2 = L
        
        # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
        N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - common.famous_f(4*D0*t/L**2) * common.famous_f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

        D_from_nmsd = 

        # fit_func = N2_func_full
        # popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
        # D0 = popt[0]
        # r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

        #, r^2={r2:.2f}
        label = rf'$L={L:.2f}\mathrm{{\mu m}}$'
        # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='grey', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)

        # linear fit to start
        # fit_end = 6
        # fit_func_2 = lambda t, D, e : 8 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D / L**2) * t**e
        # popt, pcov = scipy.optimize.curve_fit(fit_func_2, t[1:fit_end], delta_N_sq[1:fit_end])
        # ax.plot(t[1:fit_end], fit_func_2(t[1:fit_end], *popt), linestyle=':', color='gray')
        # label += rf' $D={popt[0]:.3f}, t^{{{popt[1]:.2f}}}$'

        # N2_theory = 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2)**2) # countoscope eq. 2
        # N2_theory_lowtime = 4 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t_theory) * (L + L_2) / (L * L_2)
        # ax.plot(t_theory[:LOWTIME_FIT_END], N2_theory_lowtime[:LOWTIME_FIT_END], linestyle='dashed', linewidth=1, color='black', label='sFDT (no inter.) low time' if box_size_index==0 else None)

        # p1, p2 = plateaus.calc_plateaus_for_L(sigma, phi, L)
        # ax.hlines(p1, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1, label='plateaus')
        # ax.hlines(p2, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1)

        # computed theory interactions
        # N2_theory_interactions = 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
        # ax.plot(t_theory, N2_theory_interactions, color='black', linewidth=1, label='sFDT (w/ inter.)' if box_size_index==0 else None)

        # fit to whole thing
        N2_theory = lambda t, D, N: 2 * N * (1 - common.famous_f(4*D*t/L**2)**2) # countoscope eq. 2
        fitting_points = common.exponential_integers(1, t.max())
        popt, pcov = scipy.optimize.curve_fit(N2_theory, t[fitting_points], delta_N_sq[fitting_points])
        ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='grey', linewidth=1)
        label += fr', $D_\mathrm{{fit}}={popt[0]:.3f}$'
        # ±{np.sqrt(pcov[0][0]):.3f}$'
        
        exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=label, linestyle='none', marker='o', zorder=-1)

    ax.legend(fontsize=8, loc='lower right')
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel('$t$')
    ax.set_ylabel('$\Delta N^2(t)$')
    title = f'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$'
    if not np.isnan(sigma):
        title += f', $\sigma={sigma:.3f}\mathrm{{\mu m}}$'
    title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
    ax.set_title(title)

    # N_stats = np.fromfile('../calc/Count_Data_Cpp/Exp_test_N_stats.txt', sep=' ')
    # N_stats = N_stats.reshape((-1, 5))
    # N_mean = {}
    # for i in range(N_stats.shape[0]):
    #     N_mean[int(N_stats[i, 0])] = N_stats[i, 1]
        
    ###### RESCALED PLOT ######
        
        
    # ax = rescaled_axs[mode_index]

    # for box_size_index in boxes_to_use:

    #     delta_N_sq = N2_mean[box_size_index, :]

    #     # delta_N_sq = delta_N_sq[1:] # [1:] is because the point at t=0 msd=0 plots weirdly

    #     t_over_L_sq = t/box_sizes[box_size_index]**2
    #     delta_N_sq_over_N = delta_N_sq / N_mean[box_size_index]



    #     # N2oN_func = lambda toL2, D0: 8/np.sqrt(np.pi) * np.sqrt(D0 * toL2) # countoscope eq. 3
    #     # popt, pcov = scipy.optimize.curve_fit(N2oN_func, t_over_L_sq[0:LOWTIME_FIT_END], delta_N_sq_over_N[0:LOWTIME_FIT_END])
    #     D0 = D0_from_fits[mode_index][box_size_index]

    #     # D_str = f'D={D0:.3f}±{D0_unc_from_fits[mode_index][box_size_index]:.3f} ({D0/D0_from_fits[0][box_size_index]:.2f})'
    #     f'D={D0:.3f}'

    #     ax.plot(t_over_L_sq[1:], delta_N_sq_over_N[1:], label=rf'$L={box_sizes[box_size_index]}\mathrm{{\mu m}}, {D_str}$', marker='.')
    #     # ax.plot(t_over_L_sq[0:num_timesteps//100], N2oN_func(t_over_L_sq[0:num_timesteps//100], D0))

    # t_over_L_sq = np.logspace(-2, 3)
    # # f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2
    # L = 1
    # # N2_theory = 2 * (1 - f(4 * D0 * t_over_L_sq)**theory_exponent)
    # # N2_theory_lowtime = 8 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2)
    # # ax.plot(t_over_L_sq[1:], N2_theory[1:], ':', color='black', label=f'sFDT (no inter.) D={D0:.3f}' if box_size_index==0 else None)

    # ax.legend(fontsize=9)
    # ax.semilogy()
    # ax.semilogx()
    # ax.set_xlabel('$t/L^2$')
    # ax.set_ylabel(r'$\Delta N^2(t) / \langle N \rangle$')
    # ax.set_title(f'rescaled, $\phi={phi}$, {mode} {driftremoved}')

    fig.tight_layout()
    fig.savefig(f'visualisation/figures_png/msd_{file}.png', dpi=300)
