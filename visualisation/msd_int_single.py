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

D0_from_fits     = [{}, {}]
D0_unc_from_fits = [{}, {}]

boxes_to_use = list(range(0, 8))

collapse_x = True
collapse_y = True
collapse_x = False
collapse_y = False

for file in sys.argv[1:]:
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

    data = common.load(f'intensity_counting/data/counted_{file}.npz')
    box_sizes               = data['box_sizes']
    counted_intensity_diffs = data['counted_intensity_diffs']
    avg_intensities         = data['avg_intensities']
    pixel_size              = data['pixel_size']
    time_step               = data['time_step']
    variances               = data['variances']
    particle_diameter       = data['particle_diameter']

    obs_plateaus = np.full((len(box_sizes)), np.nan)

    if file == 'sim_downsampled':
        # line up the colours
        ax.scatter([], [])
        ax.scatter([], [])

    for box_size_index in range(len(box_sizes)):
        intensity_diff = counted_intensity_diffs[box_size_index, :-201] # must be odd?
        t = np.arange(0, len(intensity_diff) * time_step, time_step)
        assert intensity_diff.shape == t.shape, f'{intensity_diff.shape} != {t.shape}'

        L = box_sizes[box_size_index]

        if collapse_y:
            collapser = avg_intensities[box_size_index]
            collapser = variances[box_size_index]
            intensity_diff                  /= collapser
            variances      [box_size_index] /= collapser
            avg_intensities[box_size_index] /= collapser
        if collapse_x:
            # t /= L**2
            # t /= L**2 + 0.160 * L + 0.215
            # xlabel = '$t/(L^2 + 0.160 * L + 0.215)$'
            t /= L**2
            xlabel = '$t/L$'
        else:
            xlabel = '$t$'

        D0 = None
        
        label = f'$L={L:.1f}\mathrm{{\mu m}}$'
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
        popt, pcov = scipy.optimize.curve_fit(fit_func_2, t[1:fit_end], intensity_diff[1:fit_end]/intensity_factor)
        # ax.plot(t[1:fit_end], fit_func_2(t[1:fit_end], *popt), linestyle=':', color='gray')
        label += rf' $D={popt[0]:.3f}, t^{{{popt[1]:.2f}}}$'

        # fit to whole thing
        t_theory = np.logspace(np.log10(t[1] / 2), np.log10(t.max()))
        N = variances[box_size_index]/intensity_factor
        N2_theory = lambda t, D: 2 * N * (1 - common.famous_f(4*D*t/L**2)**2) # countoscope eq. 2
        # ax.plot(t_theory, N2_theory(t_theory, popt[0]), color='black', linewidth=1)

        fitting_points = common.exponential_indices(t)
        popt, pcov = scipy.optimize.curve_fit(N2_theory, t[fitting_points], intensity_diff[fitting_points]/intensity_factor)
        if not collapse_x and not collapse_y:
            pass
            # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1)
        label += fr', $D_\mathrm{{fit}}={popt[0]:.3f}$'


        ax.scatter(t[1:], intensity_diff[1:]/intensity_factor, label=label, s=30, marker='.', edgecolors='none')


    # if collapse_y:
    #     predicted_plateaus /= avg_intensities[box_size_index]

    if not collapse_y:
        # ax.hlines(2*avg_intensities/np.sqrt(intensity_factor), t.max()/20, t.max(), linestyle='dashed', linewidth=1, color='black', label=r'$2\langle I \rangle / I_f{}^2$')
        ax.hlines(2*variances/intensity_factor,                t.max()/20, t.max(), linestyle='dotted', linewidth=1, color='black', label=r'$2Var(I) / I_f$')


        # ax.hlines(variances[box_size_index]/intensity_factor**2, t.max()/20, t.max(), linestyle='dotted', color='grey', label=r'$Var(I) / I_f$' if box_size_index==0 else None)


    ax.legend(loc='lower right', fontsize=6)
    # ax.legend(fontsize=5)
    # ax.semilogy()
    ax.set_ylim(0, 0.1)
    ax.set_xlim(0, 10)
    # ax.semilogx()
    ax.set_xlabel(xlabel)
    ylabel = r'$\Delta I^2(t) / \rangle I \langle I_f{}^2$' if collapse_y else '$\Delta I^2(t) / I_f{}^2$'
    ax.set_ylabel(ylabel)
    ax.set_title(fr'{file} intensity counting, $\sigma={particle_diameter}\mathrm{{\mu m}}$')

    common.save_fig(fig, f'visualisation/figures_png/msd_int_{file}.png', dpi=600)
