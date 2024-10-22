import common
import matplotlib.pyplot as plt
import numpy as np
import countoscope_theory.structure_factor
import scipy.optimize, scipy.stats

SMALL = False
SHOW_R_AXIS = False
FORCE_LOG_X_AXIS = True
RESCALE_X_AXIS_BY_DIAMETER = False
SPLIT_AND_COLOR = False
SHOW_FIT = False
SHOW_THEORY = True

def go(file, export_destination=None):
    data = common.load(f"scattering_functions/data/F_{file}.npz")
    t                 = data["t"]
    F                 = data["F"] # (num timesteps) x (num k bins)
    F_unc             = data['F_unc']
    k                 = data["k"]
    particle_diameter = data.get('particle_diameter', 1)

    # S = F[0, :, :]
    # k = k[0, :, :]
    # S_unc = F_unc[0, :, :]
    S     = F    [0, :]
    k     = k    [0, :]
    S_unc = F_unc[0, :]

    fig, ax = plt.subplots(1, 1, figsize=(3.2, 3) if SMALL else (4, 4))
    
    start_index = 40 # crops off k=0 delta fn
    start_index = 0
    end_index = 1000
    if RESCALE_X_AXIS_BY_DIAMETER:
        x = particle_diameter*k[0, :]
    else:
        x = k

    min = np.nanargmin(S)

    if SPLIT_AND_COLOR:
        pass
        # ax.errorbar(x[start_index:min], S[start_index:min], yerr=F_unc[0, start_index:min], linestyle='none', marker='o', color='tab:orange')
        # # ax.errorbar(x[start_index:min], S[start_index:min], yerr=F_unc[0, start_index:min], linestyle='none', marker='o', color='tab:green', alpha=0.5)
        # ax.errorbar(x[min:end_index],   S[min:end_index],   yerr=F_unc[0, min:end_index],   linestyle='none', marker='o', color='tab:green')
    else:
        ax.errorbar(x.flatten(), S.flatten(), yerr=S_unc.flatten(), linestyle='none', marker='o', color='tab:green')
        ax.errorbar(x.flatten(), S.flatten(), yerr=S_unc.flatten(), linestyle='none', marker='o', color='tab:green')

    if data['log'] or FORCE_LOG_X_AXIS:
        ax.semilogx()
        x_th = np.logspace(np.log10(x.min()), np.log10(x.max()))
    else:
        x_th = np.linspace(0, x.max())
    # popt, pcov = scipy.optimize.curve_fit(countoscope_theory.structure_factor.hard_spheres_2d, x[min:end_index], S[min:end_index], p0=(0.3, 2))
    # ax.plot(x_th, countoscope_theory.structure_factor.hard_spheres_2d(x_th, *popt), color='grey', label=f'fit ($\phi={common.format_val_and_unc(popt[0], np.sqrt(pcov[0, 0]), sigfigs=3)}$, $\sigma={common.format_val_and_unc(popt[1], np.sqrt(pcov[1, 1]), sigfigs=3)}$)')
    
    if SHOW_FIT:
        popt, pcov = scipy.optimize.curve_fit(lambda k, phi, s : countoscope_theory.structure_factor.hard_spheres_2d(k, phi, s), x[min:end_index], S[min:end_index], p0=(0.5, 2))
        ax.plot(x_th, countoscope_theory.structure_factor.hard_spheres_2d(x_th, *popt), color=common.FIT_COLOR, label=f'fit: $\phi={common.format_val_and_unc(popt[0], np.sqrt(pcov[0, 0]), sigfigs=3)}$, $\sigma={common.format_val_and_unc(popt[1], np.sqrt(pcov[1, 1]), sigfigs=3)}$')
    
    if SHOW_THEORY:
        if (pack_frac_given := data.get('pack_frac_given')) and (particle_diameter := data['particle_diameter']):
            ax.plot(x, countoscope_theory.structure_factor.hard_spheres_2d(x, pack_frac_given, particle_diameter))
    # ax.plot(x_th, countoscope_theory.structure_factor.hard_spheres_2d(x_th, 0.34, 3.03), color='grey', label='$\sigma=3.03$', zorder=10)


    # bins = np.linspace(k.min(), k.max(), data['num_k_bins'])
    # print(type(bins))
    # print(type(S))
    # print(type(k))
    # print(k.shape, S.shape, bins.shape)
    # F_binned = scipy.stats.binned_statistic(k.flatten(), S.flatten(), bins=bins)[0]
    # bin_mids = (bins[1:] + bins[:-1])/2
    # ax.plot(bin_mids, F_binned, label='binned', zorder=20)

    # above_2 = S>2
    # print(np.argwhere(above_2))
    # print('k', k[np.nonzero(above_2)])
    # print('k_x', data['k_x'][np.nonzero(above_2)[0]])
    # print('k_y', data['k_y'][np.nonzero(above_2)[1]])
    # print('S', S[np.nonzero(above_2)])
    # print('above 2', np.count_nonzero(above_2))
    # print('k', k[above_2])
    # print('S', S[above_2])
    # ax.plot(k[above_2], S[above_2], marker='o')


    ax.set_ylabel('$S(k)$')
    if RESCALE_X_AXIS_BY_DIAMETER:
        ax.set_xlabel('$k\sigma$')
    else:
        ax.set_xlabel('$k$')

    if SHOW_R_AXIS:
        realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k/particle_diameter, lambda r: 2*np.pi*particle_diameter/r))
        realspace_ax.set_xticks([1e1, 1e0, 5e-1])
        realspace_ax.set_xlabel(r'$r/\sigma = 2\pi/k\sigma$')



    ax.legend()
    # ax.semilogy()
    # ax.set_ylim(0.05, 10000)
    # ax.set_ylim(0, 2)

    # ax.text(0.7, 0.1, f'$\phi={file[-4:]}$', transform=ax.transAxes)

    if export_destination:
        common.save_fig(fig, export_destination, hide_metadata=True)
    common.save_fig(fig, f'scattering_functions/figures_png/S_of_k_{file}.png')



if __name__ == '__main__':
    for file in common.files_from_argv('scattering_functions/data', 'F_'):
        go(file)