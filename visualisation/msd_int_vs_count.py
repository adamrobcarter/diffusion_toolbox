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

boxes_to_use = list(range(0, 8))

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
    # print(t_theory,np.log10(t[1] / 10), np.log10(t.max()))

    for box_size_index in boxes_to_use:
    # for L in [2**e for e in range(-2, 7)]:
        L = box_sizes[box_size_index]

        delta_N_sq = N2_mean[box_size_index, :]
        # t = np.arange(0, len(delta_N_sq))[1:]/2
        # delta_N_sq = delta_N_sq # [1:] is because the point at t=0 msd=0 plots weirdly
        

        f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2
        

        L_2 = L
        
        # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
        N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc


        LOWTIME_FIT_END = 20

        fit_func = N2_func_full
        popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
        D0 = popt[0]
        r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))

        #, r^2={r2:.2f}
        D_str = f'D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}'

        exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=rf'particles $L_x={L}\mathrm{{\mu m}}, {D_str}$', linestyle='none', marker='_', markersize=3)
        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='grey', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)

        # N2_theory_interactions = 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
        # ax.plot(t_theory, N2_theory_interactions, color='black', linewidth=1, label='sFDT (w/ inter.)' if box_size_index==0 else None)



    collapse_x = False
    collapse_y = False

    data = common.load(f'intensity_counting/data/counted_{file}.npz')
    box_sizes               = data['box_sizes']
    counted_intensity_diffs = data['counted_intensity_diffs']
    avg_intensities         = data['avg_intensities']
    pixel_size              = data['pixel_size']
    time_step               = data['time_step']

    obs_plateaus = np.full((len(box_sizes)), np.nan)

    if file == 'sim_downsampled':
        # line up the colours
        ax.scatter([], [])
        ax.scatter([], [])

    for box_size_index in range(len(box_sizes)):
        intensity_diff = counted_intensity_diffs[box_size_index, :-201] # must be odd?
        t = np.arange(0, len(intensity_diff) * time_step, time_step)
        assert intensity_diff.shape == t.shape, f'{intensity_diff.shape} != {t.shape}'

        # intensity_diff = intensity_diff / avg_intensities[i]
        # t = t / box_sizes[i]**2

        # removed = intensity_diff < 1e-11
        # t = t[~removed]
        # intensity_diff = intensity_diff[~removed]

        # intensity_diff = intensity_diff / avg_intensities[box_size_index]
        # t = t / box_sizes[box_size_index]**2

        L = box_sizes[box_size_index]
        # if L > 32:
        #     continue


        if collapse_y:
            intensity_diff                  /= avg_intensities[box_size_index]
            avg_intensities[box_size_index] /= avg_intensities[box_size_index]
        if collapse_x:
            t /= L**2

        D0 = None
        
        label = f'intensity $L={L:.1f}\mathrm{{\mu m}}$'
        label += f'$= {box_sizes[box_size_index]/pixel_size:.0f}\mathrm{{px}}$'

        
        obs_plateaus[box_size_index] = np.median(intensity_diff[-100:-1])
        # nmsd_ax.hlines(obs_plateaus[box_size_index], t.max()/20, t.max(), linewidth=5, color='black', label='obs plateau' if box_size_index==0 else None)
            
        # intensity_factor = np.median(avg_intensities[box_size_index]/np.median(intensity_diff[-100:-1]))
        if file == 'pierre_sim':
            intensity_factor = 16.17**2
        else:
            intensity_factor = 10573**2
        # nmsd_ax.hlines(2 * avg_intensities[box_size_index] / intensity_factor, t.max()/20, t.max(), linestyle='dashed', color='grey', label=r'$2\langle I \rangle / I_f$' if box_size_index==0 else None)
            

        if file == 'sim' or file == 'sim_downsampled':
            print(f'sim avg/plateau (L={L:.1f}) {avg_intensities[box_size_index]/np.median(intensity_diff[-100:-1]):.0f}')

            
            print('intensity factor', intensity_factor)

            # avg_intensities[box_size_index] /= intensity_factor

            label += f' ({intensity_factor:.0f})'

            # nmsd_ax.hlines(avg_intensities[box_size_index]/intensity_factor, t.max()/20, t.max(), linestyle=':', color='grey', label='avg_int/int_factor' if box_size_index==0 else None)
            

            sigma = 0.6

        fit_end = 5
        fit_func_2 = lambda t, D, e : 8 / np.sqrt(np.pi) * avg_intensities[box_size_index]/intensity_factor * np.sqrt(D / L**2) * t**e
        popt, pcov = scipy.optimize.curve_fit(fit_func_2, t[1:fit_end], intensity_diff[1:fit_end])
        # nmsd_ax.plot(t[1:fit_end], fit_func_2(t[1:fit_end], *popt), linestyle=':', color='gray')
        label += rf' $D={popt[0]:.3f}, t^{{{popt[1]:.2f}}}$'


        ax.scatter(t[1:], intensity_diff[1:]/intensity_factor, label=label, s=10, marker='|')




        predicted_plateaus = 2 * avg_intensities / intensity_factor # one for rescaled axis
        predicted_plateaus /=  intensity_factor # one cause that's what we expect
        # predicted_plateaus /= 1.5
        
        # if collapse_y:
        #     predicted_plateaus /= avg_intensities[box_size_index]
        ax.hlines(predicted_plateaus, t.max()/20, t.max(), linestyle='dashed', color='grey', label=r'$2\langle I \rangle / I_f$' if box_size_index==0 else None)


    # print('calced', intensity_factors.find(file))

    # ax.legend(loc='lower right', fontsize=8)
    ax.legend(fontsize=5)
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel('$t$')
    ax.set_ylabel('$\Delta N^2(t)$')
    ax.set_title(f'{file} detection & counting vs intensity counting')

    fig.tight_layout()
    fig.savefig(f'visualisation/figures_png/msd_int_vs_count_{file}.png', dpi=300)
