import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_theory.nmsd, countoscope_theory.structure_factor
from .msd_single import get_plateau
import scipy.optimize

for file in common.files_from_argv('box_counting/data/', 'counted_'):
    fig, ax = plt.subplots(1, 1)
    fig2, ax2 = plt.subplots(1, 1)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean        = data['N2_mean']
    N2_std         = data['N2_std']
    phi            = data['pack_frac']
    sigma          = data['particle_diameter']
    time_step      = data['time_step']
    depth_of_field = data.get('depth_of_field')
    box_sizes      = data['box_sizes']
    N_mean         = data['N_mean']
    N_mean_std     = data['N_mean_std']
    N_var          = data['N_var']
    N_var_mod      = data['N_var_mod']
    N_var_mod_std  = data['N_var_mod_std']
    sep_sizes      = data['sep_sizes']
    time_step      = data['time_step']
    t = np.arange(0, N2_mean.shape[1]) * time_step


    num_timesteps = N2_mean.shape[1]

    obs_plateaus = np.full(box_sizes.shape, np.nan)
    obs_plateaus_unc = np.full(box_sizes.shape, np.nan)

    for i in range(len(box_sizes)):
        obs_plateaus[i], obs_plateaus_unc[i] = get_plateau('obs', nmsd=N2_mean[i, :], file=file, L=box_sizes[i], phi=phi, sigma=sigma, t=t, var=N_var[i], varmod=N_var_mod[i])

    ax.errorbar(box_sizes, obs_plateaus, yerr=0, zorder=10, label='obs', color='tab:blue')
    ax.errorbar(box_sizes, obs_plateaus, yerr=obs_plateaus_unc, alpha=0.5, color='tab:blue')
    # ax.set_xlim(sep_sizes.max()*1.1, sep_sizes.min()*1.1)
    # ax.set_xlim(16, -16)
    # ax.set_ylim(27, 35)
    # common.save_fig(fig, f'/home/acarter/presentations/countoscope_may/plateaus_justhalf_{file}.png', dpi=200)

    plat_ths = []
    for i in range(sep_sizes.size):
        plat_ths.append(
            countoscope_theory.nmsd.plateau_inter_2d(N_mean[i], box_sizes[i], lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)),
        )

    ax.plot(box_sizes, plat_ths, label='theory (inter)', color='tab:green')

    ax.set_xlabel('L')
    ax.set_ylabel('plateau')

    ax.set_title(fr'plateeaus $\phi={phi:.3f}$')

    ax.errorbar(box_sizes, N_mean*2, yerr=N_mean_std*2, alpha=1, color='tab:orange', label='2*mean')
    ax.errorbar(box_sizes, N_mean*2, yerr=0, color='tab:orange', alpha=0.5)

    ax.errorbar(box_sizes, N_var_mod*2, yerr=0,         color='tab:red', label='$2 Var_\mathrm{boxes}(N)$')
    ax.errorbar(box_sizes, N_var_mod*2, yerr=N_var_mod_std*2, color='tab:red', alpha=0.5)
    
    ax.plot(box_sizes, N_var*2, color='tab:purple', label='$2 Var(N)$')

    # p = countoscope_theory.nmsd.plateau_inter_2d(N_mean.mean(), box_sizes.mean(), lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
    # ax.hlines(p/2, box_sizes.min(), box_sizes.max(), color='tab:brown', linestyle='dashed', label='var theory')

    fitted_plats = []
    for box_size_index in range(box_sizes.size):
        L = box_sizes[box_size_index]

        def timescaleint_replacement_fitplateau(nmsd):
            N2_theory2 = lambda t, D, plateau : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
            log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
            
            fitting_points2 = common.exponential_integers(1, t.size//2)
            # p0 = (0.05, N_mean[box_size_index])
            p02 = [0.05, 2*N_var[box_size_index]]
            popt2, pcov2 = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points2], np.log(nmsd[fitting_points2]), p0=p02, maxfev=2000)
            # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
            D_from_fit2 = popt2[0]
            D_from_fit_unc2 = np.sqrt(pcov2[0][0])

            # fit_ys = N2_theory2(t_theory, *popt2)

            return D_from_fit2, D_from_fit_unc2, popt2[1]
        
        Dc, Dc_unc, plat = timescaleint_replacement_fitplateau(N2_mean[box_size_index, :])
        fitted_plats.append(plat)

    ax.plot(box_sizes, fitted_plats, color='tab:pink', label='nmsd fit')


    ax.legend(loc='center left', fontsize=8)

    ax.semilogy()
    ax.semilogx()

    if 'eleanorlong001' in file:
        ax.set_xlim(1e1, 3e2)
        ax.set_ylim(3e-1, 7e2)



    common.save_fig(fig, f'box_counting/figures_png/plateaus_{file}.png', dpi=200)

