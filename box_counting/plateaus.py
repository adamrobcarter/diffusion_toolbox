import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_theory.nmsd, countoscope_theory.structure_factor
from .msd_single import get_plateau, PLATEAU_SOURCES
import scipy.optimize
import visualisation.Ds_overlapped

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
    N_var_time     = data.get('N_var_time')
    N_var_mod_std  = data['N_var_mod_std']
    sep_sizes      = data['sep_sizes']
    time_step      = data['time_step']
    t = np.arange(0, N2_mean.shape[1]) * time_step

    D0, _, _ = visualisation.Ds_overlapped.get_D0(file)

    sources = PLATEAU_SOURCES
    sources = ['var', 'sDFT']

    for source_i, source in enumerate(sources):
        plateaus = np.full(box_sizes.shape, np.nan)
        plateaus_unc = np.full(box_sizes.shape, np.nan)

        for i in range(len(box_sizes)):
            plateaus[i], plateaus_unc[i] = get_plateau(method=source, nmsd=N2_mean[i, :], file=file,
                                L=box_sizes[i], phi=phi, sigma=sigma, t=t,
                                var=N_var[i], varmod=N_var_mod[i], D0=D0,
                                N_mean=N_mean[i], density=data.get('density', 0),
                                var_time=N_var_time, cutoff_L=box_sizes[21], cutoff_plat=2*N_var[21])
            
            # if source == 'var' and i == 21:
            #     v_10 = plateaus[i]

        ax.errorbar(box_sizes/data['particle_diameter'], plateaus, yerr=plateaus_unc, label=source, marker='.')

    # ax.errorbar(box_sizes/data['particle_diameter'], v_10 * box_sizes**2 / box_sizes[21]**2, label=r'$Var_{10\sigma} L^2/(10\sigma)^2$', marker='.')

    ax.legend()

    ax.semilogy()
    ax.semilogx()

    if 'eleanorlong001' in file:
        ax.set_xlim(1e1, 3e2)
        ax.set_ylim(3e-1, 7e2)
    if 'eleanorlong010' in file:
        ax.set_xlim(1e1, 1e2)
        ax.set_ylim(2e1, 3e3)
    if 'L640' in file or 'L320' in file:
        ax.set_xlim(1.1e1, 1.1e2)
        ax.set_ylim(2e1, 3e3)

    ax.set_ylabel('plateau')
    ax.set_xlabel('$L/\sigma$')


    common.save_fig(fig, f'box_counting/figures_png/plateaus_{file}.png', dpi=200)

