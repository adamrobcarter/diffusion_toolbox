import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common
import numpy as np

fig, ax = plt.subplots()

visualisation.Ds_overlapped_mult.go(
    (
        # dict(
        #     file='ld_hydro_nbody_0.114_L2560_t1h_1',
        #     source='f_t1',
        #     plot_index=np.index_exp[23:],
        #     # color=WALL_SHORT_COLOR,
        #     msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
        #     label='single wall, $t=1\mathrm{{s}}$'
        # ),
        # dict(
        #     file='ld_hydro_nbody_0.114_L2560',
        #     source='f_t1024',
        #     msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
        #     # color=WALL_LONG_COLOR,
        #     label='single wall, $t=1024\mathrm{{s}}$'
        # ),
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='black',
        ),
        dict(
            file='H_of_k_L2560_numkbins60_single_wall_phi0.114_maxk21.8_myRPY_numis1000',
            msd_file='ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub_2d',
            source='H',
            label='theory, singlewall hydro (numeric)',
            color='tab:red',
        ),
        # nolub
        
        # no lub
        dict(
            file='ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub',
            source='f_t1',
            msd_file='ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub_2d',
            label='no lubrication $t=1\mathrm{s}$',
            color='tab:blue',
            plot_index=np.index_exp[45:],
            # alpha = 0
        ),
        dict(
            file='ld_nbody_flat_0.114_singlewall_L2560_t8h_32s_dt50_nolub',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub_2d',
            label='no lubrication $t=1024\mathrm{s}$',
            color='tab:green',
        ),
    ),
    ax = ax,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True,
)
# ax.set_ylim(0.5, 4)
ax.set_ylim(0.8, 4)
ax.set_xlim(1e-2, 3e1)
ax.semilogy()

common.save_fig(fig, 'workflows/libm_directH.png')




print()
print()



fig, ax = plt.subplots()

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='black',
        ),
        dict(
            file='H_of_k_L2560_numkbins60_open_phi0.114_maxk21.8_myRPY_numis1000',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            source='H',
            label='theory, singlewall hydro (numeric)',
            color='tab:red',
        ),
        # nolub
        
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[19:],
            color='tab:blue',
            label='potential confinement, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t8h_32',
            source='f_t1024',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[:19],
            color='tab:green',
            label='potential confinement, $t=1024\mathrm{{s}}$'
        ),
    ),
    ax = ax,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True,
)
# ax.set_ylim(0.5, 4)
ax.set_ylim(0.8, 20)
ax.set_xlim(1e-2, 3e1)
ax.semilogy()

common.save_fig(fig, 'workflows/libm_directH_open.png')