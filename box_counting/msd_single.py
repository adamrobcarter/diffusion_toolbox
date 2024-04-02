import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
import sDFT_interactions
import sys

present_small = False
figsize = (6, 4.5)
if present_small:
    figsize = (4.5, 4)

collapse_x = True
collapse_y = True
collapse_x = False
collapse_y = False

for file in common.files_from_argv('box_counting/data/', 'counted_'):

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean        = data['N2_mean']
    N2_std         = data['N2_std']
    N_stats        = data['N_stats']
    phi            = data['pack_frac']
    sigma          = data['particle_diameter']
    time_step      = data['time_step']
    depth_of_field = data.get('depth_of_field')

    box_sizes = N_stats[:, 0]
    N_mean    = N_stats[:, 1]
    N_var     = N_stats[:, 2]
    num_of_boxes = N_stats[:, 5]

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t_all = np.arange(0, num_timesteps) * time_step

    # reduce = 1
    # t        = t      [::reduce]
    # # t_theory = np.logspace()
    # N2_mean  = N2_mean[:, ::reduce]
    # N2_std   = N2_std [:, ::reduce]


    # for box_size_index, L in enumerate(box_sizes):
    for box_size_index, L in list(enumerate(box_sizes)):
    # for L in [2**e for e in range(-2, 7)]:
        L = box_sizes[box_size_index]

        delta_N_sq     = N2_mean[box_size_index, :]
        delta_N_sq_err = N2_std [box_size_index, :]
        t = np.copy(t_all)
        t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()))

        anomalous = delta_N_sq < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        if np.any(anomalous):
            print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
            delta_N_sq     = delta_N_sq    [~anomalous]
            delta_N_sq_err = delta_N_sq_err[~anomalous]
            t              = t             [~anomalous]
        
        # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
        # N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - common.famous_f(4*D0*t/L**2) * common.famous_f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

        # fit_func = N2_func_full
        # popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
        # D0 = popt[0]
        # r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

        #, r^2={r2:.2f}
        label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
        # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        # ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='grey', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        # ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)

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
        # D0 = { # countoscope paper, table 1
        #     'alice0.02': 0.0416,
        #     'alice0.02_overlapped': 0.0416,
        #     'alice0.34': 0.0310,
        #     'alice0.66': 0.0175
        # }[file]
        # N2_theory_interactions = 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
        # ax.plot(t_theory, N2_theory_interactions, color='black', linewidth=1, label='sFDT (w/ inter.)' if box_size_index==0 else None)

        # fit to whole thing
        if depth_of_field:
            N2_theory = lambda t, D, N: common.N2_nointer_3D(t, D, N, L, L, depth_of_field)
        else:
            N2_theory = lambda t, D, N: common.N2_nointer_2D(t, D, N, L, L)
        fitting_points = np.unique(np.round(10**np.linspace(0, np.log10(t.max()))).astype('int'))
        popt, pcov = scipy.optimize.curve_fit(N2_theory, t[fitting_points], delta_N_sq[fitting_points], maxfev=2000)
        # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
        
        N2_theory_points = N2_theory(t_theory, *popt)
        
        if collapse_y:
            delta_N_sq       /= N_var[box_size_index]
            delta_N_sq_err   /= N_var[box_size_index]
            N2_theory_points /= N_var[box_size_index]
        if collapse_x:
            t /= L**2
            t_theory /= L**2
        
        
        if not collapse_y:
            ax.plot(t_theory[1:], N2_theory_points[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
        # label += fr', $D_\mathrm{{fit}}={popt[0]:.3g}\pm {np.sqrt(pcov[0][0]):.3g} \mathrm{{\mu m^2/s}}$'
        label += fr', $D_\mathrm{{fit}}={common.format_val_and_unc(popt[0], np.sqrt(pcov[0][0]), 2)} \mathrm{{\mu m^2/s}}$'
        # ±{np.sqrt(pcov[0][0]):.3f}$'
        
        if collapse_x or collapse_y:
            markersize = 2
        else:
            markersize = 5

        exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
        # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:]/np.sqrt(num_of_boxes[box_size_index]), label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
        # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
    

    ax.legend(fontsize=7, loc='lower right')
    ax.semilogy()
    ax.semilogx()
    xlabel = '$t/L^2$' if collapse_x else '$t$'
    ylabel = r'$\Delta N^2(t)/\langle N \rangle$' if collapse_y else '$\Delta N^2(t)$'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    title = file
    # title = f'Simulated colloids in RCP spheres\n$\phi={phi:.3f}$'
    if not np.isnan(phi):
        title += f', $\phi_\mathrm{{calc}}={phi:.3f}$'
    if not np.isnan(sigma):
        title += f', $\sigma={sigma:.2f}\mathrm{{\mu m}}$'
    if sigma_calced := data.get('particle_diameter_calced'):
        title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
        # print('sigma calced hidden from legend')
    ax.set_title(title)

    common.save_fig(fig, f'box_counting/figures_png/msd_{file}.png', dpi=300)
