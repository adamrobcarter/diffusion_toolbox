import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common

fig, ax = plt.subplots()

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='f_t1024',
            # plot_index=np.index_exp[23:],
            # color=WALL_SHORT_COLOR,
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            # label='single wall, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_nbody_flat_0.114_singlewall_L1280_t8h_32s_dt50',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L1280_t1h_1s_dt50_2d',
        ),
        dict(
            file='ld_nbody_flat_0.114_singlewall_L640_t8h_32s_dt125',
            source='f_t1024',
            msd_file='ld_nbody_flat_0.114_singlewall_L640_t1h_1s_dt125_2d',
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='black'
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
ax.set_ylim(0.8, 9)
ax.set_xlim(1e-2, 3e1)
ax.semilogy()

common.save_fig(fig, 'workflows/libm_var_L.png')