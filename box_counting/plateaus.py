import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_theory.nmsd, countoscope_theory.structure_factor
from .msd_single import get_plateau, PLATEAU_SOURCES, PLATEAU_SOURCE_NAMES
import scipy.optimize
import visualisation.Ds_overlapped
import tqdm

markers = {
    'sDFT': '^',
    'var': 'o',
}

def go(file, ax, sources, rescale_window=False, label_prefix='', rescale_sigma=True, colors=None,
       linestyle='-'):

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

    if colors:
        assert len(colors) == len(sources)

    for source_i, source in enumerate(sources):
        plateaus = np.full(box_sizes.shape, np.nan)
        plateaus_unc = np.full(box_sizes.shape, np.nan)

        for i in tqdm.trange(len(box_sizes)):
            if 'target' in source and i < 10:
                continue
            
            plateaus[i], plateaus_unc[i] = get_plateau(method=source, nmsd=N2_mean[i, :], file=file,
                                L=box_sizes[i], phi=phi, sigma=sigma, t=t,
                                var=N_var[i], varmod=N_var_mod[i], D0=D0,
                                N_mean=N_mean[i], density=data.get('density', 0),
                                var_time=N_var_time, cutoff_L=box_sizes[19], cutoff_plat=2*N_var[19])
            
            # if source == 'var' and i == 21:
            #     v_10 = plateaus[i]
        # ax.errorbar(box_sizes/data['particle_diameter'], plateaus, yerr=plateaus_unc, label=PLATEAU_SOUURCE_NAMES.get(source, source), marker='.')
        
        rescale_x, rescale_y = 1, 1
        if rescale_window:
            window_size = min(data['window_size_x'], data['window_size_y'])
            rescale_x *= window_size**1
            rescale_y *= window_size**2
        if rescale_sigma:
            rescale_x *= data['particle_diameter']
        
        color = None
        if colors:
            color = colors[source_i]
        
        label = label_prefix + PLATEAU_SOURCE_NAMES.get(source, source)
        ax.plot(box_sizes/rescale_x, plateaus/rescale_y, label=label, marker=markers.get(source, '.'),
                color=color, linestyle=linestyle)

    # ax.errorbar(box_sizes/data['particle_diameter'], v_10 * box_sizes**2 / box_sizes[21]**2, label=r'$Var_{10\sigma} L^2/(10\sigma)^2$', marker='.')

    ax.legend()

    ax.semilogy()
    ax.semilogx()

    if rescale_window:
        ax.set_ylabel('$\mathrm{plateau}/L_x{}^2$')
        if rescale_sigma:
            raise NotImplemented()
        else:
            ax.set_xlabel('$L/L_x$')
    else:
        ax.set_ylabel('plateau')
        ax.set_xlabel('$L/\sigma$')
        

        # if 'eleanorlong001' in file:
        #     ax.set_xlim(1e1, 3e2)
        #     ax.set_ylim(3e-1, 7e2)
        # if 'eleanorlong010' in file:
        #     ax.set_xlim(1e1, 1e2)
        #     ax.set_ylim(2e1, 3e3)
        # if 'L640' in file or 'L320' in file:
        #     ax.set_xlim(1.1e1, 1.1e2)
        #     ax.set_ylim(2e1, 3e3)


if __name__ == '__main__':
    
    for file in common.files_from_argv('box_counting/data/', 'counted_'):
        fig, ax = plt.subplots(1, 1)

        go(file, ax, ['var', 'sDFT'])
        # go(file, ax, PLATEAU_SOURCES)

        common.save_fig(fig, f'box_counting/figures_png/plateaus_{file}.png', dpi=200)