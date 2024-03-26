import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
import sys
# import sDFT_interactions
import visualisation.intensity_factors as intensity_factors

# integrate = lambda *args, **kwargs: scipy.integrate.quad(*args, **kwargs)[0]

file = sys.argv[1]

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

D0_from_fits     = [{}, {}]
D0_unc_from_fits = [{}, {}]

boxes_to_use = list(range(0, 8, 2))

collapse_x = True
collapse_y = True
collapse_x = False
collapse_y = False

for file in sys.argv[1:]:
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    data = common.load(f'box_counting/data/counted_{file}.npz')
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
    t = np.arange(0, N2_mean.shape[1] * time_step, time_step)

    t_theory = np.logspace(np.log10(t[1] / 5), np.log10(t.max()))

    for box_size_index in boxes_to_use:
        L = box_sizes[box_size_index]

        delta_N_sq = N2_mean[box_size_index, :]
    
        if collapse_y:
            collapser = N_mean[box_size_index]
            # collapser = N_var[box_size_index]
            delta_N_sq             /= collapser
            N_var [box_size_index] /= collapser
            N_mean[box_size_index] /= collapser
            
        t = np.arange(0, N2_mean.shape[1] * time_step, time_step)

        if collapse_x:
            t /= L**2

        f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2
        
        L_2 = L
        N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

        LOWTIME_FIT_END = 20

        fit_func = N2_func_full
        popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], delta_N_sq[0:LOWTIME_FIT_END])
        D0 = popt[0]
        r2 = common.r_squared(delta_N_sq[0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

        #, r^2={r2:.2f}
        D_str = f'D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}'

        print(delta_N_sq[1:])
        print(t[1:])

        # exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=rf'part $L={L}\mathrm{{\mu m}}, {D_str}$', linestyle='none', markersize=6, color='tab:blue')
        label = rf'part $L={L}\mathrm{{\mu m}}, {D_str}$'
        label = 'countoscope' if box_size_index == 0 else None
        exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=label, color='tab:blue', marker='o', linestyle='none', markersize=3)
        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        # ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='black', linestyle='dashed', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        # ax.hlines(2*N_var [box_size_index], t.min(), t.max(), color='black',                     linewidth=1, label=r'$2 Var(N)$' if box_size_index==0 else None)
        # ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)

        # N2_theory_interactions = 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
        # ax.plot(t_theory, N2_theory_interactions, color='black', linewidth=1, label='sFDT (w/ inter.)' if box_size_index==0 else None)


    data = common.load(f'intensity_counting/data/counted_{file}.npz')
    box_sizes               = data['box_sizes']
    counted_intensity_diffs = data['counted_intensity_diffs']
    avg_intensities         = data['avg_intensities']
    pixel_size              = data['pixel_size']
    time_step               = data['time_step']
    variances               = data['variances']

    obs_plateaus = np.full((len(box_sizes)), np.nan)

    if file == 'sim_downsampled':
        # line up the colours
        ax.scatter([], [])
        ax.scatter([], [])

    for box_size_index in range(0, len(box_sizes), 2):
        intensity_diff = counted_intensity_diffs[box_size_index, :-201] # must be odd?
        t = np.arange(0, len(intensity_diff) * time_step, time_step)
        assert intensity_diff.shape == t.shape, f'{intensity_diff.shape} != {t.shape}'

        L = box_sizes[box_size_index]
        # if L > 32:
        #     continue

        if collapse_y:
            collapser = avg_intensities[box_size_index]
            collapser = variances[box_size_index]
            intensity_diff                  /= collapser
            variances      [box_size_index] /= collapser
            avg_intensities[box_size_index] /= collapser
        if collapse_x:
            t /= L**2
            pass

        D0 = None
        
        label = f'int $L={L:.1f}\mathrm{{\mu m}}$'
        # label += f'$= {L/pixel_size:.0f}\mathrm{{px}}$'

        obs_plateaus[box_size_index] = np.median(intensity_diff[-100:-1])
            
        # intensity_factor = np.median(avg_intensities[box_size_index]/np.median(intensity_diff[-100:-1]))
        if file == 'pierre_sim':
            intensity_factor = 16.17**2
        else:
            intensity_factor = 10573**2
        

        if file == 'sim' or file == 'sim_downsampled':
            print(f'sim avg/plateau (L={L:.1f}) {avg_intensities[box_size_index]/np.median(intensity_diff[-100:-1]):.0f}')
            print('intensity factor', intensity_factor)
            label += f' ({intensity_factor:.0f})'

        fit_end = 5
        fit_func_2 = lambda t, D, e : 8 / np.sqrt(np.pi) * avg_intensities[box_size_index]/intensity_factor * np.sqrt(D / L**2) * t**e
        popt, pcov = scipy.optimize.curve_fit(fit_func_2, t[1:fit_end], intensity_diff[1:fit_end])
        # ax.plot(t[1:fit_end], fit_func_2(t[1:fit_end], *popt), linestyle=':', color='gray')
        # label += rf' $D={popt[0]:.3f}, t^{{{popt[1]:.2f}}}$'

        label='intensity countoscope' if box_size_index == 0 else None
        ax.scatter(t[1:], intensity_diff[1:]/intensity_factor, label=label, s=10, color='tab:orange')


        # predicted_plateaus = 2 * avg_intensities / intensity_factor # one for rescaled axis
        # predicted_plateaus /=  intensity_factor # one cause that's what we expect
        # # predicted_plateaus /= 1.5

    predicted_plateaus = 2 * avg_intensities / np.sqrt(intensity_factor)
    print(predicted_plateaus)
    # if collapse_y:
    #     predicted_plateaus /= avg_intensities[box_size_index]
    # ax.hlines(predicted_plateaus, t.max()/20, t.max(), linestyle='dashed', color='grey', label=r'$2\langle I \rangle / I_f{}^2$')


    # ax.hlines(2*variances/intensity_factor, t.max()/20, t.max(), linestyle='dotted', color='grey', label=r'$2Var(I) / I_f$')


        # ax.hlines(variances[box_size_index]/intensity_factor**2, t.max()/20, t.max(), linestyle='dotted', color='grey', label=r'$Var(I) / I_f$' if box_size_index==0 else None)


    # print('calced', intensity_factors.find(file))

    ax.legend(loc='lower right', fontsize=8)
    # ax.legend(fontsize=5)
    ax.semilogy()
    ax.semilogx()
    xlabel = '$t/L^2$' if collapse_x else '$t$'
    ax.set_xlabel(xlabel)
    ylabel = r'$\Delta N^2(t) / \rangle N \langle$' if collapse_y else '$\Delta N^2(t)$'
    ax.set_ylabel(ylabel)
    sigma = data['particle_diameter']
    ax.set_title(f'{file} countoscope vs intensity countoscope, $\sigma={sigma}\mathrm{{\mu m}}$')


    ax.text(0.4, 3e-4, '$L=0.1\mathrm{\mu m}$')
    ax.text(0.4, 1.2e-2, '$L=0.4\mathrm{\mu m}$')
    ax.text(0.4, 9e-2, '$L=1.6\mathrm{\mu m}$')
    ax.text(0.4, 6e-1, '$L=6.4\mathrm{\mu m}$')

    fig.tight_layout()
    fig.savefig(f'visualisation/figures_png/msd_int_vs_count_{file}.png', dpi=300)
