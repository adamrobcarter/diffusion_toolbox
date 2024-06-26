import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
import sDFT_interactions
import sys

# integrate = lambda *args, **kwargs: scipy.integrate.quad(*args, **kwargs)[0]

def S(k, phi, sigma):
    # sigma is disk diameter
    rho = 4 * phi / (np.pi * sigma**2)

    phi = np.pi/4 * rho * sigma**2
    J0 = lambda x: scipy.special.jv(0, x)
    J1 = lambda x: scipy.special.jv(1, x)

    prefactor = np.pi / ( 6 * ( 1 - phi)**3 * k**2 )
    line1 = -5/4 * (1 - phi)**2 * k**2 * sigma**2 * J0(k * sigma / 2)**2
    line23 = 4 * ( (phi - 20) * phi + 7) + 5/4 * (1 - phi)**2 * k**2 * sigma**2
    line23factor = J1(k * sigma / 2)**2
    line4 = 2 * (phi - 13) * (1 - phi) * k * sigma * J1(k * sigma / 2) * J0(k * sigma / 2)
    c = prefactor * (line1 + line23*line23factor + line4)
    
    S = 1 / (1 - rho * c)

    return S

file = sys.argv[1]

phi = 0.02
sigma = 2.8

D0_from_fits     = [{}, {}]
D0_unc_from_fits = [{}, {}]

LOWTIME_FIT_END = 20

boxes_to_use = list(range(0, 8))
boxes_to_use.reverse()

# driftremoved = '_driftremoved'
# driftremoved = ''

ylims = {}

# for mode_index, mode in enumerate(['boxes', 'rects']):
mode = 'boxes'
mode_index = 1

mode_option = '_long'
mode_option = ''

driftremoved_options = ['', '_driftremoved']
driftremoved_options = ['']
phis = [0.02, 0.34, 0.66]
# phis = [0.34]

for driftremoved_index, driftremoved in enumerate(driftremoved_options):

    msd_fig, msd_axs           = plt.subplots(len(phis), 1, figsize=(6, len(phis)*4.5), squeeze=False)
    rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    for phi_index, phi in enumerate(phis):
        data = common.load(f'data/counted_{mode}{mode_option}{driftremoved}_{phi}.npz')
        # data = common.load(f'data/counted_driftremoved_{phi}.npz')
        N2_mean = data['N2_mean']
        N2_std  = data['N2_std']
        N_stats = data['N_stats']
        time_step = data['time_step']

        box_sizes = N_stats[:, 0]
        N_mean    = N_stats[:, 1]
        N_var     = N_stats[:, 2]

        num_timesteps = N2_mean.shape[1]
        num_boxes     = N2_mean.shape[0]
        t = np.arange(0, num_timesteps) * time_step

        reduce = 1
        t        = t      [::reduce]
        # t_theory = np.logspace()
        N2_mean  = N2_mean[:, ::reduce]
        N2_std   = N2_std [:, ::reduce]

        t_theory = np.logspace(np.log10(t[1] / 5), np.log10(t.max()))
        # print(t_theory,np.log10(t[1] / 10), np.log10(t.max()))


        # D0 = { # countoscope paper, table 1
        #     0.02: 0.0416,
        #     0.34: 0.0310,
        #     0.66: 0.0175
        # }[phi]

        # if mode == 'strips':
        #     D0 = D0 / 2
        ax = msd_axs[phi_index][0]

        for box_size_index in boxes_to_use:
        # for L in [2**e for e in range(-2, 7)]:
            L = box_sizes[box_size_index]

            delta_N_sq = N2_mean[box_size_index, :]
            # t = np.arange(0, len(delta_N_sq))[1:]/2
            # delta_N_sq = delta_N_sq # [1:] is because the point at t=0 msd=0 plots weirdly
            
            # D0 = 0.038 # Bare diffusion coefficient in um^2/s -- short time?
            # # D0 = D0 * 2.2

            # if phi == 0.66:   
            #     D0 = D0 / 2.2 / 1.2
            #     pass

            # first_grad = ]
            f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2
            

            # ax = msd_axs[mode_index]

            if mode == 'boxes':
                L_2 = L
            elif mode == 'strips':
                L_2 = 217.6
            elif mode == 'rects':
                L_2 = 150
            
            # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
            N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc
        
            fit_func = N2_func_full
            popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
            D0 = popt[0]
            r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

            
            D0_from_fits    [mode_index][box_size_index] = D0
            D0_unc_from_fits[mode_index][box_size_index] = np.sqrt(pcov[0][0])

            #, r^2={r2:.2f}
            D_str = f'D={D0:.3f}±{D0_unc_from_fits[mode_index][box_size_index]:.3f}' if mode=='boxes' or True else f'D_s/D_b={D0/D0_from_fits[0][box_size_index]:.2f}, r^2={r2:.2f}'

            exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=rf'$L_x={L}\mathrm{{\mu m}}, {D_str}$', linestyle='none', marker='o')
            ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

            ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='grey', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
            ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)


            # N2_theory = 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2)**2) # countoscope eq. 2
            # N2_theory_lowtime = 4 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t_theory) * (L + L_2) / (L * L_2)
            # ax.plot(t_theory[:LOWTIME_FIT_END], N2_theory_lowtime[:LOWTIME_FIT_END], linestyle='dashed', linewidth=1, color='black', label='sFDT (no inter.) low time' if box_size_index==0 else None)

            # p1, p2 = plateaus.calc_plateaus_for_L(sigma, phi, L)
            # ax.hlines(p1, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1, label='plateaus')
            # ax.hlines(p2, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1)

            N2_theory_interactions = 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
            ax.plot(t_theory, N2_theory_interactions, color='black', linewidth=1, label='sFDT (w/ inter.)' if box_size_index==0 else None)

        # ax.legend(loc='lower right', fontsize=8)
        ax.legend(fontsize=8)
        ax.semilogy()
        ax.semilogx()
        ax.set_xlabel('$t$')
        ax.set_ylabel('$\Delta N^2(t)$')
        if mode == 'rects':
            title = f'rectangles, $L_y=150\mathrm{{\mu m}}$'
        else:
            title = mode
        ax.set_title(f'$\phi={phi}$, {title} {driftremoved}')

        if mode == 'strips' or mode == 'rects':
            ax.set_ylim(0.1, 10)
            pass


        if driftremoved == '':
            ylims[phi_index] = ax.get_ylim()
        elif driftremoved == '_driftremoved':
            ax.set_ylim(*ylims[phi_index])
            
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

    msd_fig     .tight_layout()
    # rescaled_fig.tight_layout()
    msd_fig     .savefig(f'figures_png/msd{driftremoved}.png',          dpi=300)
    # rescaled_fig.savefig(f'figures_png/msd{driftremoved}_rescaled.png', dpi=300)
