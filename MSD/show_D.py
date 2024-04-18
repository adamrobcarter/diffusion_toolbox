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
    ax.plot(t[1:], msd[1:]/(4*t[1:]), marker='.', markersize=3, linestyle='none', label=r'$1/4 \cdot \langle r^2 \rangle/t$')
    ax.plot(t[1:], np.gradient(msd[1:], t[1:])/4, marker='.', markersize=3, linestyle=r'none', label=r'$1/4 \cdot \mathrm{d}\langle r^2 \rangle/\mathrm{d}t$')
    # ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2)

    fitting_points = common.exponential_integers(1, t.size-1)

    func = lambda t, D: 4*D*t
    popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
    t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
    # ax.plot(t_th, func(t_th, *popt), color='black', linewidth=1, label='fit')
    ax.hlines(popt[0], t[1], t[-1], color='black')

    ax.semilogx()
    # ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    ax.set_ylabel(r'$D$')
    ax.set_xlabel('$t$ (s)')
    ax.legend()

    common.save_fig(fig, f'MSD/figures_png/D_from_msd_{file}.png')