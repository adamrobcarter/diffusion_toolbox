import matplotlib.pyplot as plt
import numpy as np
import common
import sDFT_interactions
import sys

# integrate = lambda *args, **kwargs: scipy.integrate.quad(*args, **kwargs)[0]

collapse_x = True
collapse_y = True
collapse_x = False
collapse_y = False

fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

for file_i, file in enumerate(sys.argv[1:]):

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20
    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean      = data['N2_mean']
    N2_std       = data['N2_std']
    N_stats      = data['N_stats']
    phi          = data['pack_frac']
    sigma        = data['particle_diameter']
    added_drift_x= data['added_drift_x']
    box_sizes_x = data['box_sizes_x']
    box_sizes_y = data['box_sizes_y']
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

        if box_sizes_y.size == 1:
            Ly = box_sizes_y
        else:
            Ly = box_sizes_y[box_size_index]
            
        if box_sizes_x.size == 1:
            Lx = box_sizes_x
        else:
            Lx = box_sizes_x[box_size_index]
        
        # computed theory interactions
        t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()))
        N2_theory = common.N2_nointer(t_theory, D0, N_var[box_size_index], Lx, Ly)
        
        if collapse_y:
            delta_N_sq /= N_var[box_size_index]
            N2_theory  /= N_var[box_size_index]
        if collapse_x:
            # t /= np.sqrt(Lx * Ly)
            t /= Lx * Ly
            t_theory /= Lx * Ly
            pass

        anomalous = delta_N_sq < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        if np.any(anomalous):
            print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
            delta_N_sq = delta_N_sq[~anomalous]
            t          = t         [~anomalous]


        # t = np.arange(0, len(delta_N_sq))[1:]/2
        # delta_N_sq = delta_N_sq # [1:] is because the point at t=0 msd=0 plots weirdly
        
        # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
        # N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - common.famous_f(4*D0*t/L**2) * common.famous_f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

        # fit_func = N2_func_full
        # popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
        # D0 = popt[0]
        # r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

        #, r^2={r2:.2f}
        label = rf'$L_x={Lx:.1f}\mathrm{{\mu m}}$, $L_y={Ly:.1f}\mathrm{{\mu m}}$'
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

        # fit to whole thing
        # N2_theory = lambda t, D, N: 2 * N * (1 - common.famous_f(4*D*t/L**2)**2) # countoscope eq. 2
        # fitting_points = np.unique(np.round(10**np.linspace(0, np.log10(t.max()))).astype('int'))
        # popt, pcov = scipy.optimize.curve_fit(N2_theory, t[fitting_points], delta_N_sq[fitting_points])
        # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
        # label += fr', $D_\mathrm{{fit}}={popt[0]:.3f}$'
        # ±{np.sqrt(pcov[0][0]):.3f}$'
        
        color = ['tab:blue', 'tab:orange'][file_i]
        exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=label, linestyle='none', color=color, marker='o', zorder=-1, markersize=3)
        ax.plot(t_theory, N2_theory, color='black', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

ax.legend(fontsize=7, loc='upper left')
ax.semilogy()
ax.semilogx()
ax.set_xlabel('$t$')
ax.set_ylabel('$\Delta N^2(t)$')
title = f'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$'
if not np.isnan(sigma):
    title += f', $\sigma={sigma:.3f}\mathrm{{\mu m}}$'
if sigma_calced := data.get('particle_diameter_calced'):
    title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
title += fr', $\nu_x={added_drift_x}\mathrm{{\mu m/s}}$'
ax.set_title(title)

props = dict(boxstyle='round', facecolor='grey', alpha=0.15)  # bbox features
# fig.text(1.03, 0.98, ' '.join(sys.argv), fontsize=12, verticalalignment='top', bbox=props)
command = '.'.join(sys.argv[0].split('/')[4:]).split('.py')[0] + ' ' + ' '.join(sys.argv[1:])
ax.text(-0.13, -0.15, command, transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props)


fig.tight_layout()
files = '_'.join(sys.argv[1:])
fig.savefig(f'box_counting/figures_png/msd_combined_{files}.png', dpi=300)
