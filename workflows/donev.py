import common
import matplotlib.pyplot as plt
import visualisation.Ds_overlapped_mult


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))

COLOR_FKT_SHORT = 'maroon'
COLOR_FKT_LONG = 'indianred'
COLOR_COUNTING_SHORT = 'navy'
COLOR_COUNTING_LONG = 'royalblue'

visualisation.Ds_overlapped_mult.go(
    [
        # f(k, t)
        dict(
            file='ld_hydro_nbody_0.114_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            color=COLOR_FKT_SHORT,
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='f_t1024',
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            color=COLOR_FKT_LONG
        ),
        # counting
        dict(
            file='ld_hydro_nbody_0.114_L2560_t1h_1',
            source='timescaleint_nofit_cropped_var',
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            color=COLOR_COUNTING_SHORT,
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='timescaleint_nofit_cropped_var',
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            color=COLOR_COUNTING_LONG,
        ),
    ],
    ax = ax1,
    legend_fontsize=6,
)

visualisation.Ds_overlapped_mult.go(
    [
        # f(k, t)
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            color=COLOR_FKT_SHORT,
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t8h_32',
            source='f_t1024',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            color=COLOR_FKT_LONG,
        ),
        # counting
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='timescaleint_nofit_cropped_var',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            color=COLOR_COUNTING_SHORT,
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t8h_32',
            source='timescaleint_nofit_cropped_var',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            color=COLOR_COUNTING_LONG,
        ),
    ],
    ax = ax2,
    legend_fontsize=6,
)

for ax in ax1, ax2:
    ax.semilogy()
    ax.set_ylim(0.8, 5)

common.save_fig(fig, 'workflows/figures/donev.png')

