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


fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

titles = []

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
color_index = 0
if sys.argv[1] == 'alice0.02_overlapped3' and sys.argv[2] == 'alice0.02_overlapped':
    color_index += 1
if sys.argv[1] == 'alice0.02_overlapped3' and sys.argv[2] == 'alice0.02_overlapped_neg':
    color_index += 1

for file in sys.argv[1:]:
    # ax.set_prop_cycle(None) # reset colour cycle
    color = colors[color_index]
    color_index += 1

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    boxes_to_use = list(range(0, 4)) # FIXME
    # boxes_to_use.reverse()

    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean = data['N2_mean']
    N2_std  = data['N2_std']
    N_stats = data['N_stats']
    phi     = data['pack_frac']
    sigma   = data['particle_diameter']
    time_step    = data['time_step']

    box_sizes = N_stats[:, 0]
    sep_sizes = data['sep_sizes']
    N_mean    = N_stats[:, 1]
    N_var     = N_stats[:, 2]
    num_boxes_used = N_stats[:, 5]

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
        
        L_2 = L
        
        # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
        N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

        fit_func = N2_func_full
        popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
        D0 = popt[0]
        r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

        #, r^2={r2:.2f}
        # D_str += f'±{np.sqrt(pcov[0][0]):.3f}'

        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        # ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='grey', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        # ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)


        # N2_theory = 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2)**2) # countoscope eq. 2
        # N2_theory_lowtime = 4 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t_theory) * (L + L_2) / (L * L_2)
        # ax.plot(t_theory[:LOWTIME_FIT_END], N2_theory_lowtime[:LOWTIME_FIT_END], linestyle='dashed', linewidth=1, color='black', label='sFDT (no inter.) low time' if box_size_index==0 else None)

        # p1, p2 = plateaus.calc_plateaus_for_L(sigma, phi, L)
        # ax.hlines(p1, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1, label='plateaus')
        # ax.hlines(p2, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1)
        
        D0 = 0.0416

        # N2_theory_interactions = 2 * N_var[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
        # ax.plot(t_theory, N2_theory_interactions, color='black', linestyle='dotted', linewidth=1, zorder=3, label='sFDT (w/ inter.)' if box_size_index==0 else None)

        label = rf'$L = {L}\mathrm{{\mu m}}$'
        # label += f', D={D0:.3f}'
        label += f', $sep = {sep_sizes[box_size_index]:.1f}\mathrm{{\mu m}}$'
        label += f', $n = {num_boxes_used[box_size_index]:.0f}$'

        ax.plot(t[1:], delta_N_sq[1:], label=label, color=color, linestyle='none', marker='.', markersize=3, markeredgecolor='none')
        
    titles.append(f'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$, $\sigma={sigma}$')

    # ax.legend(loc='lower right', fontsize=8)
legend = ax.legend(fontsize=6)
for handle in legend.legend_handles:
    handle.set_markersize(6.0)
ax.semilogy()
ax.semilogx()
ax.set_xlabel('$t$')
ax.set_ylabel('$\Delta N^2(t)$')
ax.set_title(', '.join(titles))
ax.set_title(titles[0])


fig.tight_layout()
names = '_'.join(sys.argv[1:])
fig.savefig(f'visualisation/figures_png/msd_combined_{names}.png', dpi=300)
