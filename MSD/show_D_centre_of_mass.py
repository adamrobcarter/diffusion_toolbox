import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize
import matplotlib.cm

for file in common.files_from_argv('MSD/data', 'msd_centre_of_mass_'):
    fig, ax = plt.subplots(1, 1)

    data = common.load(f'MSD/data/msd_centre_of_mass_{file}.npz')

    for group_index, groupsize in enumerate(groupsizes := data['groupsizes']):
        msd = data['msds'][group_index, :]
        msd_unc = data['msd_uncs'][group_index, :]


        t = np.arange(0, msd.size) * data['time_step']

        # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
        color =  matplotlib.cm.afmhot(np.interp(group_index, (0, len(groupsizes)), (0.2, 0.75)))
        ax.plot(t[1:], msd[1:]/(4*t[1:]), marker='.', color=color, markersize=3, linestyle='none', label=rf'g={groupsize} $1/4 \cdot \langle r^2 \rangle/t$')
        # ax.plot(t[1:], np.gradient(msd[1:], t[1:])/4, color=color, marker='.', markersize=3, linestyle=r'none', label=rf'g={groupsize} $1/4 \cdot \mathrm{d}\langle r^2 \rangle/\mathrm{d}t$')
        # ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2)

        fitting_points = common.exponential_integers(1, t.size-1)

        func = lambda t, D: 4*D*t
        popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
        t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
        # ax.plot(t_th, func(t_th, *popt), color='black', linewidth=1, label='fit')
        ax.hlines(popt[0], t[1], t[-1], color='black')

    # ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$D$')
    ax.set_xlabel('$t$ (s)')
    ax.legend()
    # ax.semilogx()
    # ax.semilogy()

    common.save_fig(fig, f'MSD/figures_png/D_from_msd_centre_of_mass_{file}.png')