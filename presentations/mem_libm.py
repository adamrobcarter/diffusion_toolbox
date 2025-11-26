import common
import particle_detection.show_movie
import matplotlib.pyplot as plt
import visualisation.Ds_overlapped_mult
import numpy as np
import matplotlib.ticker

path = '/home/acarter/presentations/mem_libm/figures'

# full size video
# particle_detection.show_movie.go(
#     'ld_hydro_nbody_open_0.114_L2560',
#     outfile = f'{path}/movie_fullsize',
#     # crop = 100,
#     every_nth_frame=10,
#     highlights = False,
#     output_type='frames',
#     figsize_mult = 2,
#     dpi = 150,
#     particle_color='black',
#     annotation_color='black'
# )

ylim = (0.8, 6)
xlim = (0.58e-1, 0.2e2)

fig, ax = plt.subplots(figsize=(3, 3.1))
ax.semilogy()
visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file='sim_nohydro_011_L1280',
            source='f_first_first',
            msd_file='sim_nohydro_011_L1280',
            # plot_index=np.index_exp[19:],
            # color=OPEN_SHORT_COLOR,
            # label='potential confinement, $t=1\mathrm{{s}}$'
            label = 'simulation - short time',
            marker = 'o',
            plot_index=np.index_exp[19:],
        ),
        dict(
            file='sim_nohydro_011_L1280_longer',
            source='f_first_first',
            msd_file='sim_nohydro_011_L1280',
            # plot_index=np.index_exp[19:],
            # color=OPEN_SHORT_COLOR,
            # label='potential confinement, $t=1\mathrm{{s}}$'
            label = 'simulation - long time',
            marker = 'o',
            plot_index=np.index_exp[:22],

        ),
        dict(
            file='sim_nohydro_011_L1280',
            source='D0Sk_theory',
            label = 'theory',
            color = 'black'
        )
    ],
    ax = ax,
    plot_against_k = True,
    show_twin_k_axis = False,
    allow_rescale_x=True,
    legend_fontsize=9,
)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)

# Major ticks at powers of 10
ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
# Minor ticks exist but are unlabeled
ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.yaxis.labelpad = -10

common.save_fig(fig, f'{path}/theory_nohydro.pdf', hide_metadata=True)

fig, ax = plt.subplots(figsize=(3, 3.1))
ax.semilogy()
visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[17:],
            # color=OPEN_SHORT_COLOR,
            # label='potential confinement, $t=1\mathrm{{s}}$'
            label = 'sim. - short time'
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t8h_32',
            source='f_t1024',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[1:20],
            # color=OPEN_LONG_COLOR,
            # label='potential confinement, $t=1024\mathrm{{s}}$'
            label = 'sim. - long time'
        ),
        dict(
            file='H_of_k_L2560_numkbins60_open_phi0.114_maxk10.0_theory_zpos0.00a_zwidth0.0a_rcutoff1.0_stepmult1.00_order3',
            source='H_theory',
            label = 'theory',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            color = 'black'
        ),
        # dict(
        #     file='H_of_k_L2560_numkbins60_open_phi0.114_maxk10.0_theory_zpos0.00a_zwidth1.0a_rcutoff1.0_stepmult1.00_order3',
        #     source='H_theory',
        #     label = 'theory 1.0a'
        # ),
    ],
    ax = ax,
    plot_against_k = True,
    show_twin_k_axis = False,
    allow_rescale_x=True,
    legend_fontsize=9,
)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)


# Major ticks at powers of 10
ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
# Minor ticks exist but are unlabeled
ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.yaxis.labelpad = -10

# common.add_exponential_index_indicator(ax, -1, (3e-1, 5), 'k', (1e-2, 1e0))
common.save_fig(fig, f'{path}/theory_open.pdf', hide_metadata=True)


lub_msd_file = 'ld_hydro_nbody_0.114_L2560_t1h_1_2d'
fig, ax = plt.subplots(figsize=(3, 3.1))
ax.semilogy()
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)


# Major ticks at powers of 10
ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
# Minor ticks exist but are unlabeled
ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.yaxis.labelpad = -10




visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file='ld_hydro_nbody_0.114_L2560_t1h_1',
            source='f_t1',
            # plot_index=np.index_exp[23:],
            plot_index=np.index_exp[15:],
            # color='cornflowerblue',
            msd_file=lub_msd_file,
            # label='lubrication $t=1\mathrm{s}$',
            label = 'simulation - short time',
            # alpha=0,
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='f_t1024',
            msd_file=lub_msd_file,
            # color='tab:blue',
            # label='lubrication $t=1024\mathrm{s}$',
            # label='single wall, $t=1024\mathrm{{s}}$'
            label = 'simulation - long time',
        ),
    ],
    ax = ax,
    plot_against_k = True,
    show_twin_k_axis = False,
    allow_rescale_x=True,
    legend_fontsize=9,
    # allow_rescale_y=False
)
common.save_fig(fig, f'{path}/sim_wall.pdf', hide_metadata=True)


visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file='H_of_k_L2560_numkbins100_single_wall_lubrication_phi0.114_maxk10.0_theory_zpos1.06a_zwidth0.0a_rcutoff1.0_stepmult1.00_order3',
            source='H_theory',
            msd_file=lub_msd_file,
            label='theory',
            # color='navy',
            linestyle='dashed',
            color = 'black'
        ),
    ],
    ax = ax,
    plot_against_k = True,
    show_twin_k_axis = False,
    allow_rescale_x=True,
    legend_fontsize=9,
    # allow_rescale_y=False
)

# common.add_exponential_index_indicator(ax, -1, (3e-1, 5), 'k', (1e-2, 1e0))
common.save_fig(fig, f'{path}/theory_wall.pdf', hide_metadata=True)





import MSD.show
fig, ax = plt.subplots(figsize=(3.5, 3))
MSD.show.go(
    'eleanorlong066_div8',
    ax,
    SHOW_SHORT_FIT=False,
    crop_end=-1000
)
common.save_fig(fig, f'{path}/msd_shortlong.png', hide_metadata=True, dpi=300)

# fig, ax = plt.subplots(figsize=(3.5, 3))
# MSD.show.go(
#     'eleanorlong066_div8',
#     ax,
#     SHOW_SHORT_FIT=False,
#     crop_end=-1000
# )
# common.save_fig(fig, f'{path}/msd_shortlong.pdf', hide_metadata=True)