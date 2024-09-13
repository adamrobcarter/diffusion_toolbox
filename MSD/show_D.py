import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

for file in common.files_from_argv('MSD/data', 'msd_'):
    data = common.load(f'MSD/data/msd_{file}.npz')
    msd = data['msd']
    msd_unc = data['msd_unc']

    t = np.arange(0, msd.size)
    n = (msd.size-1) / t
    msd_unc = msd_unc / np.sqrt(n)
    # ^^ this is based on thinking about how many independent measurements there are

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3))

    t = np.arange(0, msd.size) * data['time_step']

    END = -200

    # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
    # ax.plot(t[1:], msd[1:]/(4*t[1:]), marker='.', markersize=3, linestyle='none', label=r'$1/4 \cdot \langle r^2 \rangle/t$')
    D = np.gradient(msd, t)/4
    D_unc = D * msd_unc/msd
    ax.plot(t[1:END], D[1:END], marker='.', markersize=3, linestyle=r'none', label=r'$1/4 \cdot \mathrm{d}\langle r^2 \rangle/\mathrm{d}t$')
    # ax.fill_between(t[1:END], D[1:END]-D_unc[1:END], D[1:END]+D_unc[1:END], alpha=0.2)

    D_smooth = moving_average(D, 100)
    # ax.plot(t[100+50:END-49], D_smooth[100:END], marker='.', markersize=3, linestyle=r'none')
    

    fitting_points = common.exponential_integers(1, t.size-1)

    func = lambda t, D: 4*D*t
    popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
    t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
    # ax.plot(t_th, func(t_th, *popt), color='black', linewidth=1, label='fit')
    # ax.hlines(popt[0], t[1], t[-1], color='black')

    ax.semilogx()
    # ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
    # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

    # ax.set_ylim(0.01, 0.04)
    ax.set_ylim(np.quantile(D, 0.02), np.quantile(D, 0.98))

    ax.set_ylabel(r'$D$')
    ax.set_xlabel('$\Delta t$ (s)')
    # ax.legend()

    # common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/D_from_msd_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/D_from_msd_{file}.png')