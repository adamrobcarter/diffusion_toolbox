import visualisation.Ds_overlapped
import matplotlib.pyplot as plt
import common
import numpy as np

fig, ax = plt.subplots(1, 1)

Ds = []
D_uncs = []
Lxs = []

D_std_msds = []

for file in (files := common.files_from_argv('visualisation/data', 'Ds_from_MSD_first_')):
    data = common.load(f'visualisation/data/Ds_from_MSD_first_centre_of_mass_entire_{file}.npz')
    # N = 
    D_CoM = data['Ds'][0]
    D_CoM_unc = data['D_uncs'][0]
    window_size_x = data['window_size_x']
    window_size_y = data['window_size_y']
    assert window_size_x == window_size_y
    N_density = data['density'] * window_size_x**2
    N = data['N_particles']
    if window_size_x > 100:
        assert np.isclose(N_density, N, rtol=0.01)
    Ds.append(D_CoM*N)
    D_uncs.append(D_CoM_unc*N)
    Lxs.append(window_size_x)

    D_std_msd, _, _ = visualisation.Ds_overlapped.get_D0(file)
    D_std_msds.append(D_std_msd)

ax.errorbar(Lxs, D_std_msds, linestyle='none', marker='o', label='standard MSD')
ax.errorbar(Lxs, Ds,         linestyle='none', marker='o', label='N * CoM MSD')
ax.loglog()

ax.set_ylim(np.mean(D_std_msd)-0.001, np.mean(D_std_msd)+0.001)

ax.set_xlabel('$L_x$')
ax.set_ylabel('D')
ax.legend()

filestring = '_'.join(files)[:50]
common.save_fig(fig, f'visualisation/figures_png/Ds_msd_window_{filestring}.png')

