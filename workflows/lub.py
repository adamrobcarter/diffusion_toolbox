import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common
import numpy as np
import visualisation.Ds_overlapped

fig, ax= plt.subplots(figsize=(5, 4))

nolub_msd_file = 'ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub_2d'
lub_msd_file = 'ld_hydro_nbody_0.114_L2560_t1h_1_2d'

lub_D0, _, _ = visualisation.Ds_overlapped.get_D0(lub_msd_file)
nolub_D0, _, _ = visualisation.Ds_overlapped.get_D0(nolub_msd_file)

rescale_y = False

visualisation.Ds_overlapped_mult.go(
    (
        # no lub
        dict( # D(k) short time
            file='ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub',
            source='f_t1',
            msd_file=nolub_msd_file,
            label='no lubrication $t=1\mathrm{s}$',
            color='orange',
            plot_index=np.index_exp[45:],
            # alpha = 0
        ),
        dict( # D(k) long time
            file='ld_nbody_flat_0.114_singlewall_L2560_t8h_32s_dt50_nolub',
            source='f_t1024',
            msd_file=nolub_msd_file,
            label='no lubrication $t=1024\mathrm{s}$',
            color='tab:orange',
        ),
        dict(
            file='H_of_k_L2560_numkbins100_single_wall_phi0.114_maxk10.0_theory_zpos1.15a_zwidth0.0a_rcutoff1.0_stepmult1.00_order3',
            source='H_theory',
            msd_file=nolub_msd_file,
            label='no lubrication theory',
            color='sienna',
        ),
        dict( # msd
            file='ld_nbody_flat_0.114_singlewall_L2560_t1h_1s_dt50_nolub',
            source='MSD_first',
            label='no lubrication MSD $t=1\mathrm{s}$',
            color='orange',
            msd_file=nolub_msd_file,
        ),

        # lub
        dict(
            file='ld_hydro_nbody_0.114_L2560_t1h_1',
            source='f_t1',
            plot_index=np.index_exp[23:],
            color='cornflowerblue',
            msd_file=lub_msd_file,
            label='lubrication $t=1\mathrm{s}$',
            # alpha=0,
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='f_t1024',
            msd_file=lub_msd_file,
            color='tab:blue',
            label='lubrication $t=1024\mathrm{s}$',
            # label='single wall, $t=1024\mathrm{{s}}$'
        ),
        dict(
            file='H_of_k_L2560_numkbins100_single_wall_lubrication_phi0.114_maxk10.0_theory_zpos1.06a_zwidth0.0a_rcutoff1.0_stepmult1.00_order3',
            source='H_theory',
            msd_file=nolub_msd_file,
            label='lubrication theory',
            color='navy',
        ),
        dict( # msd
            file='ld_hydro_nbody_0.114_L2560_t1h_1',
            source='MSD_first',
            color='cornflowerblue',
            msd_file=lub_msd_file,
            label='lubrication MSD $t=1\mathrm{s}$',
        ),

        # theory
        # dict(
        #     file='ld_hydro_nbody_0.114_L2560',
        #     source='D0Sk_theory',
        #     msd_file=lub_msd_file,
        #     color='black',
        #     label='theory, no hydrodynamics',
        # ),
    ),
    ax = ax,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True,
    allow_rescale_y=rescale_y
)
# ax1.set_ylim(0.5, 4)
if rescale_y:
    ax.set_ylim(0.8, 5)

    ax.text(3e-2, 1.3, f'$D_0={lub_D0  :.3g}\\mathrm{{\\mu m^2/s}}$', color='tab:blue',   fontsize=12)
    ax.text(3e-2, 1.1, f'$D_0={nolub_D0:.3g}\\mathrm{{\\mu m^2/s}}$', color='tab:orange', fontsize=12)
else:
    ax.set_ylim(0.03, 0.3)

ax.set_xlim(1e-2, 3e1)
ax.semilogy()

# ax.set_title('$L=640\mathrm{\mu m}$')
common.save_fig(fig, 'workflows/figures/libm_lub_nolub.png')