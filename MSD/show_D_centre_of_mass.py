import matplotlib.pyplot as plt
import common
import numpy as np
import scipy.optimize

SHOW_REGULAR_MSD = False

for file in common.files_from_argv('MSD/data', 'msd_centre_of_mass_'):
    fig, ax = plt.subplots(1, 1)
    data = common.load(f'MSD/data/msd_centre_of_mass_{file}.npz')
    t = data['t']
    print('t', t)
    print(data['groupsizes'])

    for group_index, groupsize in enumerate(groupsizes := data['groupsizes']):
        print('doing', groupsize)
        msd               = data['msds'][group_index, :]
        msd_unc           = data['msd_uncs'][group_index, :]
        density           = data['density']
        particle_diameter = data['particle_diameter']


        color = common.colormap(group_index, 0, len(groupsizes))

        f = msd*groupsize
        y = np.gradient(f, t) / 4

        # ax.errorbar(t, msd*groupsize, msd_unc*groupsize, linestyle='none', marker='none', color=color)
        label = fr'$N={groupsize}$, $L\approx{np.sqrt(groupsize/density)/particle_diameter:.2g}\sigma$'
        ax.errorbar(t, y, marker='.', markersize=8, linestyle='none', color=color, label=label)
        # ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2)


    # show from the regular MSD calculation
    if SHOW_REGULAR_MSD:
        data_msd = common.load(f'MSD/data/msd_{file}.npz')
        msd = data_msd['msd']
        t = np.arange(0, msd.size) * data_msd['time_step']
        D = np.gradient(msd, t)/4
        ax.plot(t, D, marker='.', markersize=3, linestyle='none', label='regular MSD', color='tab:blue')




    ax.loglog()

    ax.set_ylabel(r'$N\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
    ax.set_xlabel('$t$ (s)')
    ax.legend()

    # common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/msd_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'MSD/figures_png/D_centre_of_mass_{file}.png')
    # np.savez(f'visualisation/data/Ds_from_MSD_{file}',
    #          Ds=[popt[0]], D_uncs=[np.sqrt(pcov)[0][0]], labels=[''])
    # np.savez(f'visualisation/data/Ds_from_MSD_short_{file}',
    #          Ds=[popt_short[0]], D_uncs=[np.sqrt(pcov_short)[0][0]], labels=[''])
    # np.savez(f'visualisation/data/Ds_from_MSD_long_{file}',
    #          Ds=[popt_long[0]], D_uncs=[np.sqrt(pcov_long)[0][0]], labels=[''])


# import matplotlib.pyplot as plt
# import common
# import numpy as np
# import scipy.optimize
# import matplotlib.cm

# for file in common.files_from_argv('MSD/data', 'msd_centre_of_mass_'):
#     fig, ax = plt.subplots(1, 1)

#     data = common.load(f'MSD/data/msd_centre_of_mass_{file}.npz')

#     for group_index, groupsize in enumerate(groupsizes := data['groupsizes']):
#         msd = data['msds'][group_index, :]
#         msd_unc = data['msd_uncs'][group_index, :]


#         t = np.arange(0, msd.size) * data['time_step']

#         # ax.errorbar(t[1:], msd[1:], msd_unc[1:], linestyle='none', marker='none', color='lightskyblue')
#         color =  matplotlib.cm.afmhot(np.interp(group_index, (0, len(groupsizes)), (0.2, 0.75)))
#         ax.plot(t[1:], msd[1:]/(4*t[1:]), marker='.', color=color, markersize=3, linestyle='none', label=rf'g={groupsize} $1/4 \cdot \langle r^2 \rangle/t$')
#         # ax.plot(t[1:], np.gradient(msd[1:], t[1:])/4, color=color, marker='.', markersize=3, linestyle=r'none', label=rf'g={groupsize} $1/4 \cdot \mathrm{d}\langle r^2 \rangle/\mathrm{d}t$')
#         # ax.fill_between(t[1:], msd[1:]-msd_unc[1:], msd[1:]+msd_unc[1:], alpha=0.2)

#         fitting_points = common.exponential_indices(t)

#         func = lambda t, D: 4*D*t
#         popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], msd[fitting_points])
#         t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
#         # ax.plot(t_th, func(t_th, *popt), color='black', linewidth=1, label='fit')
#         ax.hlines(popt[0], t[1], t[-1], color='black')

#     # ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
#     # ax.set_xlim(t[1]*0.8, t[-1]/0.8)

#     ax.set_ylabel(r'$D$')
#     ax.set_xlabel('$t$ (s)')
#     ax.legend()
#     # ax.semilogx()
#     # ax.semilogy()

#     common.save_fig(fig, f'MSD/figures_png/D_from_msd_centre_of_mass_{file}.png')