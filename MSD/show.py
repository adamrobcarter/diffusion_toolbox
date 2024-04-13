import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize

for file in common.files_from_argv('MSD/data', 'msd_'):
    data = np.load(f'MSD/data/msd_{file}.npz')
    msd = data['msd']
    msd_unc = data['msd_unc']

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3))

    t = np.arange(0, msd.size) * data['time_step']

    # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
    ax.plot(t[1:], msd[1:], marker='.', markersize=8, linestyle='none', color='tab:blue', label='observations')
    # ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2)

    fitting_points = common.exponential_integers(1, t.size-1)

    func = lambda t, D: 4*D*t
    popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
    t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
    ax.plot(t_th, func(t_th, *popt), color='black', linewidth=1, label='fit')

    ax.loglog()
    ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$t$ (s)')
    ax.legend()

    common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/msd_{file}.png', hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/msd_{file}.png')
    np.savez(f'visualisation/data/Ds_from_MSD_{file}',
             Ds=[popt[0]], D_uncs=[np.sqrt(pcov)[0][0]], labels=[''])