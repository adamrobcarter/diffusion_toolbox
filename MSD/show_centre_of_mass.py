import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize
import matplotlib.cm

for file in common.files_from_argv('MSD/data', 'msd_centre_of_mass_'):
    fig, ax = plt.subplots(1, 1)
    data = common.load(f'MSD/data/msd_centre_of_mass_{file}.npz')
    for group_index, groupsize in enumerate(groupsizes := data['groupsizes']):
        print('doing', groupsize)
        msd = data['msds'][group_index, :]
        msd_unc = data['msd_uncs'][group_index, :]


        t = np.arange(0, msd.size) * data['time_step']

        color =  matplotlib.cm.afmhot(np.interp(group_index, (0, len(groupsizes)), (0.2, 0.75)))

        print(common.nanfrac(msd))

        # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
        ax.plot(t[1:], msd[1:], marker='.', markersize=8, linestyle='none', color=color, label=f'g={groupsize}')
        # ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2)

        fitting_points = common.exponential_integers(1, t.size-1)
        func = lambda t, D, c: 4*D*t**c
        log_func = lambda t, D, c: np.log10(func(t, D, c))
        if common.nanfrac(msd[fitting_points]):
            print('skipping')
            continue
        popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
        # popt, pcov = scipy.optimize.curve_fit(log_func, np.log10(t[fitting_points]), np.log10(msd[fitting_points]))
        t_th = np.logspace(np.log10(t[fitting_points[0]]), np.log10(t[fitting_points[-1]]))
        ax.plot(t_th, func(t_th, *popt), color='black', linewidth=1, label=f'fit $r^2 = 4 * {popt[0]:.3f} * t^{{{popt[1]:.2f}}}$')

        # fitting_points_short = common.exponential_integers(1, t.size//100)
        # func_short = lambda t, D: 4*D*t
        # popt_short, pcov_short = scipy.optimize.curve_fit(func_short, t[fitting_points_short], msd[fitting_points_short])
        # t_th_short = np.logspace(np.log10(t[fitting_points_short[0]]), np.log10(t[fitting_points_short[-1]]))
        # ax.plot(t_th_short, func_short(t_th_short, *popt_short), color='black', linewidth=1, label='fit')

        # fitting_points_long = common.exponential_integers(t.size//10, t.size-1)
        # func_long = lambda t, D, a: 4*D*t + a
        # popt_long, pcov_long = scipy.optimize.curve_fit(func_long, t[fitting_points_long], msd[fitting_points_long])
        # t_th_long = np.logspace(np.log10(t[fitting_points_long[1]]), np.log10(t[fitting_points_long[-1]]))
        # ax.plot(t_th_long, func_long(t_th_long, *popt_long), color='black', linewidth=1, label='fit')
        # break

    ax.loglog()
    # ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$t$ (s)')
    ax.legend()

    # common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/msd_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/msd_centre_of_mass_{file}.png')
    # np.savez(f'visualisation/data/Ds_from_MSD_{file}',
    #          Ds=[popt[0]], D_uncs=[np.sqrt(pcov)[0][0]], labels=[''])
    # np.savez(f'visualisation/data/Ds_from_MSD_short_{file}',
    #          Ds=[popt_short[0]], D_uncs=[np.sqrt(pcov_short)[0][0]], labels=[''])
    # np.savez(f'visualisation/data/Ds_from_MSD_long_{file}',
    #          Ds=[popt_long[0]], D_uncs=[np.sqrt(pcov_long)[0][0]], labels=[''])