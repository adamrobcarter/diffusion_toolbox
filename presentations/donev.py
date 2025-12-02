import common
import matplotlib.pyplot as plt
import visualisation.Ds_overlapped_mult
import matplotlib.image
import matplotlib.offsetbox
import numpy as np

path = '/home/acarter/presentations/donev/figures'


WALL_LONG_COLOR = 'royalblue'
WALL_SHORT_COLOR = 'cornflowerblue'
OPEN_LONG_COLOR = 'indianred'
OPEN_SHORT_COLOR = 'lightcoral'
THEORY_COLOR = 'indigo'
WALLTRAP_LONG_COLOR = 'seagreen'
WALLTRAP_SHORT_COLOR = 'mediumseagreen'

FIG_DOUBLE_WIDTH = 7
FIG_SINGLE_WIDTH = FIG_DOUBLE_WIDTH / 2
FIG_SINGLE_HEIGHT = 3.2

def show_png_on_axis(ax, file, coords, size, color='black', use_data_coords=True):
    image = matplotlib.image.imread(file)
    imagebox = matplotlib.offsetbox.OffsetImage(image, zoom=size)   

    ab = matplotlib.offsetbox.AnnotationBbox(
        imagebox,
        coords,        # position in DATA units
        xycoords="data" if use_data_coords else 'axes fraction', # <-- use data coordinates
        # frameon=False,
    )
    
    frame = ab.patch
    frame.set_edgecolor(color)     # frame/border color
    # frame.set_linewidth(2)         # thickness

    ax.add_artist(ab)

####### main result ##########
fig, ax = plt.subplots(figsize=(5, 4))
visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_0.114_L2560_t1h_1',
            source='f_t1',
            plot_index=np.index_exp[23:],
            color=WALL_SHORT_COLOR,
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            label='single wall, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='f_t1024',
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d',
            color=WALL_LONG_COLOR,
            label='single wall, $t=1024\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[19:],
            color=OPEN_SHORT_COLOR,
            label='potential confinement, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t8h_32',
            source='f_t1024',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[:19],
            color=OPEN_LONG_COLOR,
            label='potential confinement, $t=1024\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color=THEORY_COLOR,
        ),
    ),
    ax = ax,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True
)
# ax.set_ylim(0.5, 4)
ax.set_ylim(0.8, 8.2)
ax.set_xlim(6e-2, 3e1)
# ax.set_xlim(2.1e-3, 1e1)
ax.semilogy()
common.add_exponential_index_indicator(ax, -1, (2.3e-1, 5e0), 'k', x_limits=(3e-2, 3.5e-1), )
# show_png_on_axis(ax, 'presentations/libmobility/side_view_2.png', (12e0, 3.5), 0.04, color=WALL_LONG_COLOR)
# show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (12e0, 6), 0.04, color=OPEN_LONG_COLOR)
show_png_on_axis(ax, 'presentations/libmobility/side_view_2.png', (0.88, 0.42), 0.05, color=WALL_LONG_COLOR, use_data_coords=False)
show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (0.88, 0.62), 0.05, color=OPEN_LONG_COLOR, use_data_coords=False)

ax.text(0.6, 0.65, '$\phi=0.11$', transform=ax.transAxes)

common.save_fig(fig, f'{path}/result.pdf', hide_metadata=True, dpi=300)
common.save_fig(fig, f'{path}/result.png', hide_metadata=True, dpi=300)


########## literature D(k)s #########
# fig, ax = plt.subplots(figsize=(3, 3))

# # this data is from Segre & Pusey 1996 fig 1b
# # QR = np.array([0.862,1.031,1.196,1.369,1.551,1.724,1.88,2.219,2.392,2.524,2.647,2.755,2.837,2.953,3.027,3.068,3.109,3.159,3.249,3.299,3.365,3.414,3.513,3.546,3.604,3.703,3.843,4.033,4.181,4.602,5.097,5.757,6.054,6.309,7.212])
# # DsQ_ov_D0 = 1/np.array([0.698,0.726,0.754,0.811,0.868,0.94,0.954,1.181,1.395,1.637,1.865,2.135,2.662,2.947,3.359,3.587,3.972,4.413,5.409,5.936,6.292,6.648,6.875,6.406,5.893,5.552,4.783,4.456,3.587,3.288,3.786,4.384,4.199,4.641,4.399])
# # # L_ov_R = 2*np.pi/QR
# # ax.scatter(QR, DsQ_ov_D0, label='Segrè and Pusey 1996, $\phi=0.47$')

# # this is Banchio et al 2018 fig 6 (top)
# q = np.array([0.004,0.005,0.005,0.006,0.006,0.007,0.007,0.008,0.008,0.009,0.009,0.01,0.011,0.011,0.012,0.012,0.013,0.013,0.013,0.014,0.014,0.015,0.015,0.015,0.015,0.016,0.016,0.016,0.016,0.017,0.017,0.017,0.017,0.017,0.018,0.018,0.018,0.018,0.018,0.019,0.019,0.019,0.019,0.019,0.019,0.02,0.02,0.021,0.021,0.021,0.022,0.022,0.023,0.023,0.024,0.024,0.025,0.025,0.026,0.026,0.027,0.027,0.028,0.028])
# Dsq_ov_d0 = 1/np.array([0.088,0.088,0.105,0.105,0.123,0.123,0.123,0.14,0.14,0.149,0.167,0.175,0.211,0.263,0.281,0.281,0.368,0.439,0.491,0.579,0.711,0.886,1.053,1.149,1.342,1.561,1.746,2,2.193,2.439,2.649,2.719,2.807,2.754,2.623,2.456,2.333,2.114,1.921,1.825,1.719,1.623,1.509,1.421,1.351,1.246,1.193,1.132,1.088,1.061,0.991,0.974,0.965,0.965,0.982,1.009,1.009,1.044,1.105,1.158,1.219,1.246,1.386,1.43])
# r = 136
# # L = 2*np.pi/q / r
# ax.scatter(q*r, Dsq_ov_d0, label='Banchio et al 2018, $\phi=0.15$')

# # Riese et al 2000: might actually show Dc plateau
# q        = np.array([0.006,0.008,0.011,0.014,0.016,0.019,0.021,0.023,0.025,0.027,0.028,0.032,0.034,0.036,0.038,0.04, 0.041,0.043,0.047])
# # ^ invers nm
# Dq_ov_D0 = 1/np.array([0.23, 0.201,0.208,0.245,0.245,0.253,0.29, 0.401,0.401,0.52, 0.699,1.19, 1.97, 2.55, 3.086,3.249,3.123,2.788,1.866])
# # this is short time, what are the others?
# r = 55 # nm
# # phi = 0.089
# print(len(q), len(Dq_ov_D0))
# ax.scatter(q*r, Dq_ov_D0, label='Riese et al 2000, $\phi=0.09$')

# ax.loglog()
# ax.set_xlabel('$kr$')
# ax.set_ylabel('$D_c(k)/D_0$')
# ax.legend(fontsize=10)

# common.save_fig(fig, f'{path}/literature_3d.pdf', hide_metadata=True)


####### trap and wall ######
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_DOUBLE_WIDTH, FIG_SINGLE_HEIGHT))
visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_0.114_singlewall_ztrap_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_0.114_singlewall_ztrap_L2560_t1h_1_2d',
            plot_index=np.index_exp[47:],
            color=WALLTRAP_SHORT_COLOR,
            label='wall & potential confinement, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_singlewall_ztrap_L2560_t8h_32',
            source='f_t1024',
            msd_file='ld_hydro_nbody_0.114_singlewall_ztrap_L2560_t1h_1_2d',
            # plot_index=np.index_exp[:19],
            color=WALLTRAP_LONG_COLOR,
            label='wall & potential confinement, $t=1024\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color=THEORY_COLOR,
        ),
    ),
    ax = ax1,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True
)

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_0.114_singlewall_ztrap_w0p5a_h8a_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_0.114_singlewall_ztrap_w0p5a_h8a_L2560_t1h_1_2d',
            plot_index=np.index_exp[47:],
            color=WALLTRAP_SHORT_COLOR,
            label='wall & potential confinement, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_singlewall_ztrap_w0p5a_h8a_L2560_t1h_1',
            source='f_t64',
            msd_file='ld_hydro_nbody_0.114_singlewall_ztrap_w0p5a_h8a_L2560_t1h_1_2d',
            plot_index=np.index_exp[:46],
            color=WALLTRAP_LONG_COLOR,
            label='wall & potential confinement, $t=64\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color=THEORY_COLOR,
        ),
    ),
    ax = ax2,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True
)


# ax.set_ylim(0.5, 4)
ax1.set_ylim(0.8, 8.2)
ax1.set_xlim(6e-2, 3e1)
# ax.set_xlim(2.1e-3, 1e1)
ax1.semilogy()


# ax.set_ylim(0.5, 4)
ax2.set_ylim(0.8, 8.2)
ax2.set_xlim(6e-2, 3e1)
# ax.set_xlim(2.1e-3, 1e1)
ax2.semilogy()
# common.add_exponential_index_indicator(ax, -1, (2.3e-1, 5e0), 'k', x_limits=(3e-2, 3.5e-1), )
# show_png_on_axis(ax, 'presentations/libmobility/side_view_2.png', (12e0, 3.5), 0.04, color=WALL_LONG_COLOR)
# show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (12e0, 6), 0.04, color=OPEN_LONG_COLOR)
show_png_on_axis(ax1, '/home/acarter/presentations/donev/figures/side_view-potential_and_wall-annotated.drawio.png', (0.8, 0.6), size=0.05, color=WALLTRAP_LONG_COLOR, use_data_coords=False)
show_png_on_axis(ax2, '/home/acarter/presentations/donev/figures/side_view-potential_and_wall-higher-annotated.drawio.png', (0.8, 0.6), size=0.05, color=WALLTRAP_LONG_COLOR, use_data_coords=False)
# show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (0.88, 0.6), 0.04, color=OPEN_LONG_COLOR, use_data_coords=False)

common.save_fig(fig, f'{path}/potentialandwall.pdf', hide_metadata=True, dpi=300)

########### periodic non periodic ##########
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_DOUBLE_WIDTH, FIG_SINGLE_HEIGHT))
ax1.semilogy()
ax2.semilogy()

COLOR_WINDOWED = 'tab:blue'
COLOR_NONWINDOWED = 'tab:pink'

visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file = 'sim_nohydro_002_L320',
            source = 'f_first_first',
            label = 'periodic boundary',
            marker='o',
        ),
        dict(
            file = 'sim_nohydro_002_L640_crop320',
            source = 'f_first_first',
            label = 'non periodic boundary',
            marker='o',
            color='tab:pink',
        ),
    ],
    ax=ax1,
    plot_against_k=True,
    allow_rescale_x=True,
    show_twin_k_axis=False,
    legend_fontsize=8,
    # markers=MARKER_FKT,
    # colors=[['olivedrab'], ['darkgreen']]
)

visualisation.Ds_overlapped_mult.go(
    # ['eleanorlong001', 'eleanorlong010'],
    [
        dict(
            file='eleanorlong010',
            source='D0Sk_theory',
            label='theory (no hydro)',
            color=THEORY_COLOR,
        ),
        dict(
            file='eleanorlong010',
            source='f_first_first',
            label = 'experiment',
            color=COLOR_NONWINDOWED,
            marker='o',
        ),
        dict(
            file='eleanorlong010_bhwindow',
            source='f_first_first',
            label= 'exp., BH window',
            color=COLOR_WINDOWED,
            marker='o',
        ),
    ],
    # labels=['theory', 'exp.', 'exp., BH window'],
    ax=ax2,
    plot_against_k=True,
    allow_rescale_x=True,
    show_twin_k_axis=False,
    # **Ds_ov_mult_props,
    # markers=['none', MARKER_FKT, MARKER_FKT],
    # colors=[COLOR_PHI011_FKT_THEORY, COLOR_PHI011_FKT, 'tab:blue'],
    legend_fontsize=8,
)
ax2.set_ylim(0.87, 4)
ax1.set_ylim(0.87, 40)
# ax_b.yaxis.set_major_locator(ticks_0p5)
DS_OVERLAPPED_XLIM_K = (0.26, 45)
ax1.set_xlim(0.1, 45)
ax2.set_xlim(0.26, 45)

show_png_on_axis(ax2, 'presentations/window_viz.png', (0.75, 0.4), size=0.04, use_data_coords=False, color=COLOR_WINDOWED)

ax1.text(0.7, 0.7, '$\phi=0.02$', transform=ax1.transAxes)
ax2.text(0.7, 0.6, '$\phi=0.11$', transform=ax2.transAxes)

common.add_exponential_index_indicator(ax1, -2, (0.3, 10), 'k', x_limits=(0.01, 0.4), )


common.save_fig(fig, f'{path}/periodic_effects.pdf', hide_metadata=True, dpi=300)



########### two walls ##########
COLOR_TWOWALLS_SHORT = 'mediumturquoise'
COLOR_TWOWALLS_LONG = 'lightseagreen'
fig, ax = plt.subplots(figsize=(FIG_SINGLE_WIDTH, FIG_SINGLE_HEIGHT))
visualisation.Ds_overlapped_mult.go(
    # ['eleanorlong001', 'eleanorlong010'],
    [
        dict(
            file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
            plot_index=np.index_exp[58:],
            color=COLOR_TWOWALLS_SHORT,
            label='two walls, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1',
            source='f_t64',
            msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
            plot_index=np.index_exp[:58],
            color=COLOR_TWOWALLS_LONG,
            label='two walls, $t=64\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_dpstokes_0.114_twowalls_L1280_t8h_32s',
            source='f_t64',
            msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
            plot_index=np.index_exp[:58],
            # color=COLOR_TWOWALLS_LONG,
            label='two walls, $t=64\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_dpstokes_0.114_twowalls_L1280_t8h_32s',
            source='f_t256',
            msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
            plot_index=np.index_exp[:58],
            # color=COLOR_TWOWALLS_LONG,
            label='two walls, $t=256\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_dpstokes_0.114_twowalls_L1280_t8h_32s',
            source='f_t1024',
            msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
            plot_index=np.index_exp[:58],
            # color=COLOR_TWOWALLS_LONG,
            label='two walls, $t=1024\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color=THEORY_COLOR,
        ),
    ],
    # labels=['theory', 'exp.', 'exp., BH window'],
    ax=ax,
    plot_against_k=True,
    allow_rescale_x=True,
    show_twin_k_axis=False,
    # **Ds_ov_mult_props,
    # markers=['none', MARKER_FKT, MARKER_FKT],
    # colors=[COLOR_PHI011_FKT_THEORY, COLOR_PHI011_FKT, 'tab:blue'],
    legend_fontsize=8,
)
# ax.set_ylim(0.8, 6)
# ax.set_xlim(6e-2, 3e1)
ax.set_ylim(0.8, 6)
ax.set_xlim(1e-2, 1e2)
ax.semilogy()

show_png_on_axis(ax, 'presentations/side_view-two_walls-1p5sigma.drawio.png', (0.75, 0.46), size=0.05, use_data_coords=False, color=COLOR_TWOWALLS_LONG)

ax.text(0.7, 0.63, '$\phi=0.11$', transform=ax.transAxes)


common.save_fig(fig, f'{path}/twowalls.pdf', hide_metadata=True, dpi=300)