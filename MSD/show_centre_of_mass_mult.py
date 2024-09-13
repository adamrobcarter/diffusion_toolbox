import matplotlib.pyplot as plt
import common
import numpy as np

fig, ax = plt.subplots(1, 1)

TARGET_L = 74*2.82

files = common.files_from_argv('MSD/data', 'msd_centre_of_mass_')
for file_i, file in enumerate(files):
    data = common.load(f'MSD/data/msd_centre_of_mass_{file}.npz')
    t                 = data['t']
    density           = data['density']
    particle_diameter = data['particle_diameter']
    groupsizes        = data['groupsizes']
    msds              = data['msds']
    msd_uncs          = data['msd_uncs']

    group_index = 30
    # groupsize = groupsizes[group_index]
    # msd        = msds[group_index, :]
    # msd_unc           = msd_uncs[group_index, :]


    color = ['tab:blue', 'tab:orange', 'tab:green'][file_i]

    # print('L', np.sqrt(groupsize/density), np.sqrt(groupsize/density)/particle_diameter)
    # # ax.errorbar(t, msd*groupsize, msd_unc*groupsize, linestyle='none', marker='none', color=color)
    # label = file + fr' $N={groupsize}$, $L\approx{np.sqrt(groupsize/density)/particle_diameter:.2g}\sigma$'
    # print(t.shape, msd)
    # ax.errorbar(t, msd*groupsize, yerr=msd_unc*groupsize, marker='.', markersize=8, linestyle=':', color=color, label=label)

    for group_index, groupsize in enumerate(groupsizes := data['groupsizes']):
        L = np.sqrt(groupsize/data['density'])
        if L < TARGET_L:
            continue
        else:
            print('going in ', L, group_index)

            print('L', np.sqrt(groupsize/density), np.sqrt(groupsize/density)/particle_diameter)
            # ax.errorbar(t, msd*groupsize, msd_unc*groupsize, linestyle='none', marker='none', color=color)
            label = file + fr' $N={groupsize}$, $L\approx{np.sqrt(groupsize/density)/particle_diameter:.2g}\sigma$'
            ax.errorbar(t, msds[group_index, :]*groupsize, marker='.', markersize=8, linestyle=':', color=color, label=label)
            break

ax.loglog()
# ax.set_ylim(msd[1:].min()*0.6, msd.max()/0.8)
# ax.set_xlim(t[1]*0.8, t[-1]/0.8)

ax.set_ylabel(r'$N\langle r(t)^2 \rangle$ ($\mathrm{\mu m}$)')
ax.set_xlabel('$\Delta t$ (s)')
ax.legend()

filename = '_'.join(files)
common.save_fig(fig, f'MSD/figures_png/msd_centre_of_mass_mult_{filename}.png')
