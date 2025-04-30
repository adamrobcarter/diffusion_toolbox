import matplotlib.pyplot as plt
import numpy as np
import common
import countoscope_theory.nmsd, countoscope_theory.structure_factor
from .msd_single import get_plateau, PLATEAU_SOURCES, PLATEAU_SOURCE_NAMES
import scipy.optimize
import visualisation.Ds_overlapped
import tqdm
from .D_of_L import T_integrand_func, timescaleintegral_nofit
import countoscope_theory
import joblib

memory = joblib.Memory('/data2/acarter/toolbox/box_counting/cache', verbose=0)

markers = {
    'sDFT': 'None',
    'var': 'o',
    'fit': 'x',
}
linestyles = {
    'sDFT': '-',
    'var': 'none', # None?
}
SOURCES = ['var', 'sDFT', 'fit']
SOURCES = ['var', 'sDFT']

CROP_OFF_L_BELOW = 10

@memory.cache
def get_Ds_over_grid(try_plateaus, L, t, nmsd, nmsd_std):
    Ds = np.full_like(try_plateaus, np.nan)

    for try_plateau_i, try_plateau in enumerate(tqdm.tqdm(try_plateaus)):
        T_integrand, T_integrand_unc = T_integrand_func(nmsd, nmsd_std, try_plateau, 0)
        T_integrand_max = T_integrand + T_integrand_unc
        T_integrand_min = T_integrand - T_integrand_unc
        D_of_L_nofit2, D_of_L_min_nofit2, D_of_L_max_nofit2 = timescaleintegral_nofit(t, L, T_integrand, T_integrand_min, T_integrand_max)
        Ds[try_plateau_i] = D_of_L_nofit2
    return Ds

def go(file, ax, sources, rescale_window=False, label_prefix='', rescale_sigma=True, colors=None):

    if 'zoom' in file:
        CROP_OFF_L_BELOW = 0

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
    N_var_losecorr = data.get('N_var_losecorr')
    N_var_mod_std  = data['N_var_mod_std']
    sep_sizes      = data['sep_sizes']
    time_step      = data['time_step']
    t = np.arange(0, N2_mean.shape[1]) * time_step

    D0, _, _ = visualisation.Ds_overlapped.get_D0(file)

    if colors:
        assert len(colors) == len(sources)

    all_plateaus = []

    for source_i, source in enumerate(sources):

        # if source == 'sDFT':
        #     # this is a weird hack so sDFT can have more Ls than box sizes in the counted data
        #     box_sizes_used = np.logspace(np.log10(box_sizes.min()), np.log10(box_sizes.max()), num=100)

        plateaus = np.full(box_sizes.shape, np.nan)
        plateaus_unc = np.full(box_sizes.shape, np.nan)

        for box_size_index in tqdm.trange(len(box_sizes)):
            if 'target' in source and box_size_index < 10:
                continue

            L = box_sizes[box_size_index]

            def get():
                if source == 'sDFT':
                    # this is a weird hack so sDFT can have more Ls than box sizes in the counted data
                    N_mean2 = data['density'] * L**2
                    return countoscope_theory.nmsd.plateau_inter_2d(N_mean2, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)), 0

                elif source == 'fit':
                    assert 'nohydro' in file

                    try_plateaus = np.logspace(np.log10(N_var[box_size_index]*1), np.log10(N_var[box_size_index]*10), num=500)

                    target_D = countoscope_theory.timescaleint.D_of_L(L, D0, phi, sigma)

                    Ds = get_Ds_over_grid(try_plateaus, L, t, N2_mean[box_size_index, :], N2_std[box_size_index, :])

                    closest_index = np.argmin(np.abs(Ds - target_D))
                    closest_plateau = try_plateaus[closest_index]
                    print('D dff', target_D - Ds[closest_index], closest_index)
                    return closest_plateau, 0


                else:
                    return get_plateau(
                        method=source, file=file, data=data, box_size_index=box_size_index,
                        phi=phi, sigma=sigma, t=t, D0=D0,
                    )
        
            plateaus[box_size_index], plateaus_unc[box_size_index] = get()
            
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

        print(source, plateaus_unc/plateaus)
        
        if source == 'sDFT':
            label = PLATEAU_SOURCE_NAMES.get(source, source)
            ax.plot(box_sizes/rescale_x, plateaus/rescale_y, label=label, marker=markers.get(source, '.'),
                    color=color, linestyle=linestyles.get(source, 'none'))
        else:
            label = label_prefix + PLATEAU_SOURCE_NAMES.get(source, source)
            ax.errorbar(box_sizes/rescale_x, plateaus/rescale_y, plateaus_unc/rescale_y, label=label, marker=markers.get(source, '.'),
                    color=color, linestyle=linestyles.get(source, 'none'))

        all_plateaus = np.concatenate([all_plateaus, plateaus])

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
        
    if CROP_OFF_L_BELOW:
        first_index = np.argmax(box_sizes > CROP_OFF_L_BELOW)
        ax.set_xlim(box_sizes[first_index]/rescale_x, box_sizes[-1]/rescale_x*1.1)
        ax.set_ylim(plateaus[first_index]/rescale_y, max(all_plateaus)/rescale_y*1.1)

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

        go(file, ax, SOURCES)
        # go(file, ax, PLATEAU_SOURCES)

        common.save_fig(fig, f'box_counting/figures_png/plateaus_{file}.png', dpi=200)