import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize
import matplotlib.cm

def go(file, SHOW_ERRORBARS=False, SHOW_FIT=True, SHOW_SHORT_FIT=False, SHOW_LONG_FIT=False, export_destination=None):
    data = common.load(f'MSD/data/msd_{file}.npz')
    msd = data['msd']
    msd_unc = data['msd_unc']
    
    t_indexes = np.arange(0, msd.size)
    n = (msd.size-1) / t_indexes
    msd_unc = msd_unc / np.sqrt(n)
    # ^^ this is based on thinking about how many independent measurements there are

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3))

    t = np.arange(0, msd.size) * data['time_step']

    # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
    color='tab:blue'
    ax.plot(t[1:], msd[1:], marker='.', markersize=8, linestyle='none', color=color, label='observations')
    if SHOW_ERRORBARS:
        ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2, color=color)
    
    ax.loglog()
    ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$\langle r(\Delta t)^2 \rangle$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$\Delta t$ (s)')
    
    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/msd_nofit_{file}.pdf', hide_metadata=True)

    fitting_points = common.exponential_integers(1, t.size-1)
    func = lambda t, D: 4*D*t
    popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
    t_th = np.logspace(np.log10(t[fitting_points[0]]), np.log10(t[fitting_points[-1]]))
    if SHOW_FIT:
        ax.plot(t_th, func(t_th, *popt), color='white', linewidth=1, label='fit')

    fitting_points_short = common.exponential_integers(1, 10)
    func_short = lambda t, D: 4*D*t
    popt_short, pcov_short = scipy.optimize.curve_fit(func_short, t[fitting_points_short], msd[fitting_points_short])
    t_th_short = np.logspace(np.log10(t[fitting_points_short[0]]), np.log10(t[fitting_points_short[-1]]))
    if SHOW_SHORT_FIT:
        ax.plot(t_th_short, func_short(t_th_short, *popt_short), color='white', linewidth=1)

    fitting_points_long = common.exponential_integers(t.size//10, t.size-1)
    func_long = lambda t, D, a: 4*D*t + a
    popt_long, pcov_long = scipy.optimize.curve_fit(func_long, t[fitting_points_long], msd[fitting_points_long])
    t_th_long = np.logspace(np.log10(t[fitting_points_long[1]]), np.log10(t[fitting_points_long[-1]]))
    if SHOW_LONG_FIT:
        ax.plot(t_th_long, func_long(t_th_long, *popt_long), color='white', linewidth=1)
    
    first_point_D = msd[1] / (2 * 2 * t[1])
    first_point_D_unc = msd_unc[1] / (2 * 2 * t[1])

    print(f'full  fit: D={popt[0]:.4f}')
    print(f'short fit: D={popt_short[0]:.4f}')
    print(f'long  fit: D={popt_long[0]:.4f}')
    print(f'first p  : D={first_point_D:.4f}')

    ax.legend()

    filename = f'msd_{file}'
    if SHOW_FIT:
        filename += '_fit'
    if export_destination:
        common.save_fig(fig, f'{export_destination}/{filename}.pdf', hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/{filename}.png')
    common.save_data(f'visualisation/data/Ds_from_MSD_{file}',
             Ds=[popt[0]], D_uncs=[np.sqrt(pcov)[0][0]], labels=[''])
    common.save_data(f'visualisation/data/Ds_from_MSD_short_{file}',
             Ds=[popt_short[0]], D_uncs=[np.sqrt(pcov_short)[0][0]], labels=[''])
    common.save_data(f'visualisation/data/Ds_from_MSD_long_{file}',
             Ds=[popt_long[0]], D_uncs=[np.sqrt(pcov_long)[0][0]], labels=[''])
    common.save_data(f'visualisation/data/Ds_from_MSD_first_{file}',
             Ds=[first_point_D], D_uncs=[first_point_D_unc], labels=[''])

if __name__ == '__main__':
    for file in common.files_from_argv('MSD/data', 'msd_'):
        go(file)
        