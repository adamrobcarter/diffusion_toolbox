import common
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)

D = 0
D_TIMES_N = 1
FUNCTION = D_TIMES_N

for file in (files := common.files_from_argv('MSD/data', 'msd_centre_of_mass_onepoint')):
    data = common.load(f'MSD/data/msd_centre_of_mass_onepoint_{file}.npz')
    msds       = data['msds']
    msd_uncs   = data['msd_uncs']
    groupsizes = data['groupsizes']
    time_step  = data['time_step']

    print(msds)


    Ds = msds / (2 * 2 * time_step)
    D_uncs = msd_uncs / (2 * 2 * time_step)

    x = np.sqrt(groupsizes / data['density'])
    if FUNCTION == D:
        y = Ds
        yerr = D_uncs
        ylabel = 'D'
    elif FUNCTION == D_TIMES_N:
        y = Ds * groupsizes
        yerr = D_uncs * groupsizes
        ylabel = 'D * N'

    ax.errorbar(x, y, yerr=yerr, label=file, linestyle='none', marker='o')

    common.save_data(f'visualisation/data/Ds_from_MSD_centre_of_mass_onepoint_{file}',
        Ds=y, D_uncs=yerr, Ns=groupsizes, density=data.get('density'),
        pixel_size=data['pixel_size'], window_size_x=data['window_size_x'], window_size_y=data['window_size_y'],
    )

filenames = '_'.join(files)

ax.set_xlabel('$\sqrt{N/\phi}$')
ax.set_ylabel(ylabel)
ax.legend()
ax.semilogy()
ax.semilogx()
common.save_fig(fig, f'MSD/figures_png/D_of_N_{filenames}.png')