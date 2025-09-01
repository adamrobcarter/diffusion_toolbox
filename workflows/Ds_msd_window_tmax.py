import visualisation.Ds_overlapped
import matplotlib.pyplot as plt
import common
import numpy as np

fig, ax = plt.subplots(1, 1)


interaction_type = 'hydro'
# interaction_type = 'nohydro'
# interaction_type = 'nointer'
interaction_names = {
    'nointer': 'no interactions',
    'nohydro': 'steric interactions',
    'hydro': 'with hydro'
}

tmaxs = {
    'nointer': ['1e4', '1e5', '1e6', '1e7'],
    'nohydro': ['1e4', '1e5', '1e6'],
    'hydro': ['1e4'],
}[interaction_type]

phi = 0.1

for tmax_i, tmax in enumerate(tmaxs):

    Ds = []
    D_uncs = []
    Lxs = []
    D_std_msds = []
    
    for file in (files := common.files_from_filenames(
        'visualisation/data',
        'Ds_from_MSD_first_centre_of_mass_entire_',
        [f'sim_{interaction_type}_{phi}_L*_t{tmax}_unwrap']
    )):
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

    color = common.colormap(tmax_i, 0, 2)

    ax.errorbar(Lxs, D_std_msds, linestyle='none', color=color, marker='x', label=rf'$t_\mathrm{{max}}={tmax}$ standard MSD')
    ax.errorbar(Lxs, Ds,         linestyle='none', color=color, marker='o', label=rf'$t_\mathrm{{max}}={tmax}$ N * CoM MSD')
ax.loglog()

# ax.set_ylim(0.035, 0.045)

ax.set_xlabel('simulation box $L_x$')
ax.set_ylabel('D')
ax.legend(fontsize=8, loc='upper left')

ax.set_title(interaction_names[interaction_type] + f', $\phi={phi}$')

filestring = interaction_type
common.save_fig(fig, f'visualisation/figures_png/Ds_msd_window_tmax_{filestring}.png')

