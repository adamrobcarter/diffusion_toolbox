import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize
import matplotlib.cm

SHOW_NOLOG_SHORTTIME = False
SHOW_ERRORBARS = True

def go(file, show_errorbars=False, SHOW_FIT=False, SHOW_SHORT_FIT=True, SHOW_LONG_FIT=False, export_destination=None):
    data = common.load(f'MSD/data/msd_{file}.npz')
    msd = data['msd']
    msd_unc = data['msd_unc']
    t = data['t']

    print('nan', common.nanfrac(msd))
    
    t_indexes = np.arange(0, msd.size)
    n = (msd.size-1) / t_indexes
    msd_unc = msd_unc / np.sqrt(n)
    # ^^ this is based on thinking about how many independent measurements there are

    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))

    # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
    color='tab:blue'
    ax.plot(t[1:], msd[1:], marker='.', markersize=8, linestyle='none', color=color, label='observations')
    if show_errorbars:
        ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2, color=color)
    
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)
    # ax.set_xlim(0, 20)
    # ax.set_ylim(0, 0.04)
    if SHOW_NOLOG_SHORTTIME:
        END = min(5000, t.size-1)
        ax.set_xlim(0, t[END])
        ax.set_ylim(0, msd[END])
    else:
        ax.loglog()
        ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)

    ax.set_ylabel(r'$\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$t$ (s)')

    print(f'<x>({t[1]}) = {np.sqrt(msd[1])/0.288:.3g} * 0.288um, <x>({t[32]}) = {np.sqrt(msd[32])/0.288} * 0.288um')
    
    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/msd_nofit_{file}.pdf', hide_metadata=True)
    
    fits = fit_msd(t, msd, msd_unc)

    print('first D=' + common.format_val_and_unc(fits['first']['D'], fits['first']['D_unc'], sigfigs=3))

    if t.size > 100:
        if SHOW_FIT:
            ax.plot(fits['full']['t'], fits['full']['MSD'], color=common.FIT_COLOR, linewidth=1, label='fit')
        print('fit D=' + common.format_val_and_unc(fits['full']['D'], fits['full']['D_unc'], sigfigs=3))

        if SHOW_SHORT_FIT:
            ax.plot(fits['short']['t'], fits['short']['MSD'], color=common.FIT_COLOR, linewidth=1, label='short fit')
        print('fit short D=' + common.format_val_and_unc(fits['short']['D'], fits['short']['D_unc'], sigfigs=3))

        if SHOW_LONG_FIT:
            ax.plot(fits['long']['t'], fits['long']['MSD'], color=common.FIT_COLOR, linewidth=1, label='long fit')
        print('fit long D=' + common.format_val_and_unc(fits['long']['D'], fits['long']['D_unc'], sigfigs=3))

    ax.legend()
    ax.grid(alpha=0.3)

    filename = f'msd_{file}'
    if SHOW_FIT:
        filename += '_fit'
    if export_destination:
        common.save_fig(fig, export_destination, hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/{filename}.png')


    metadata = dict(
        particle_diameter=data.get('particle_diameter'),
        pack_frac_given  =data.get('pack_frac_given'),
        pack_frac        =data.get('pack_frac'),
        pixel_size       =data.get('pixel_size'),
        window_size_x    =data.get('window_size_x'),
        window_size_y    =data.get('window_size_y')
    )

    if t.size > 100:
        common.save_data(f'visualisation/data/Ds_from_MSD_{file}',
            Ds=[fits['full']['D']], D_uncs=[fits['full']['D_unc']], labels=[''],
            **metadata
        )
        common.save_data(f'visualisation/data/Ds_from_MSD_short_{file}',
            Ds=[fits['short']['D']], D_uncs=[fits['short']['D_unc']], labels=[''],
            **metadata
        )
        common.save_data(f'visualisation/data/Ds_from_MSD_long_{file}',
            Ds=[fits['long']['D']], D_uncs=[fits['long']['D_unc']], labels=[''],
        **metadata
        )
    common.save_data(f'visualisation/data/Ds_from_MSD_first_{file}',
        Ds=[fits['first']['D']], D_uncs=[fits['first']['D_unc']], labels=[''],
        **metadata
    )

def fit_msd(t, msd, msd_unc):
    ret = {}

    if t.size > 100:
        ####### full fit #######
        fitting_points = common.exponential_indices(t)
        func = lambda t, D: 4*D*t
        popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])

        t_th = np.logspace(np.log10(t[fitting_points[0]]), np.log10(t[fitting_points[-1]]))
        fit = func(t_th, *popt)
        ret['full'] = {
            'D': popt[0],
            'D_unc': np.sqrt(pcov[0, 0]),
            't': t_th,
            'MSD': fit,
        }
        # print(f'full  fit: D={popt[0]:.4f}')

        ######### short fit ##########
        fitting_points_short = range(1, min(10, t.size-1))
        func_short = lambda t, D: 4*D*t
        popt_short, pcov_short = scipy.optimize.curve_fit(func_short, t[fitting_points_short], msd[fitting_points_short])
        
        t_th_short = np.logspace(np.log10(t[fitting_points_short[0]]), np.log10(t[fitting_points_short[-1]]))
        fit_short = func_short(t_th_short, *popt_short)
        ret['short'] = {
            'D': popt_short[0],
            'D_unc': np.sqrt(pcov_short[0, 0]),
            't': t_th_short,
            'MSD': fit_short,
        }
        # print(f'short fit: D={popt_short[0]:.4f}')

        ######### long fit ###########
        fitting_points_long = common.exponential_integers(t.size//10, t.size-1)
        func_long = lambda t, D, a: 4*D*t + a
        popt_long, pcov_long = scipy.optimize.curve_fit(func_long, t[fitting_points_long], msd[fitting_points_long])
        
        t_th_long = np.logspace(np.log10(t[fitting_points_long[1]]), np.log10(t[fitting_points_long[-1]]))
        fit_long = func_long(t_th_long, *popt_long)
        ret['long'] = {
            'D': popt_long[0],
            'D_unc': np.sqrt(pcov_long[0, 0]),
            't': t_th_long,
            'MSD': fit_long,
        }
        # print(f'long  fit: D={popt_long[0]:.4f}')

    assert msd[1] > 0
    first_point_D = msd[1] / (2 * 2 * t[1])
    first_point_D_unc = msd_unc[1] / (2 * 2 * t[1])
    ret['first'] = {
        'D': first_point_D,
        'D_unc': first_point_D_unc
    }
    # print(f'first p  : D={first_point_D:.4f}')

    return ret

if __name__ == '__main__':
    for file in common.files_from_argv('MSD/data', 'msd_'):
        go(file, show_errorbars=SHOW_ERRORBARS)
        