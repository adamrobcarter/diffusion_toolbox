import common
import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import presentations.presentations
import numpy as np

path = '/home/acarter/presentations/csi2/figures'


########### periodic non periodic ##########
fig, ax1 = plt.subplots(figsize=(3, 3))
ax1.semilogy()

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

# visualisation.Ds_overlapped_mult.go(
#     # ['eleanorlong001', 'eleanorlong010'],
#     [
#         dict(
#             file='eleanorlong010',
#             source='D0Sk_theory',
#             label='theory (no hydro)',
#             color=THEORY_COLOR,
#         ),
#         dict(
#             file='eleanorlong010',
#             source='f_first_first',
#             label = 'experiment',
#             color=COLOR_NONWINDOWED,
#             marker='o',
#         ),
#         dict(
#             file='eleanorlong010_bhwindow',
#             source='f_first_first',
#             label= 'exp., BH window',
#             color=COLOR_WINDOWED,
#             marker='o',
#         ),
#     ],
#     # labels=['theory', 'exp.', 'exp., BH window'],
#     ax=ax2,
#     plot_against_k=True,
#     allow_rescale_x=True,
#     show_twin_k_axis=False,
#     # **Ds_ov_mult_props,
#     # markers=['none', MARKER_FKT, MARKER_FKT],
#     # colors=[COLOR_PHI011_FKT_THEORY, COLOR_PHI011_FKT, 'tab:blue'],
#     legend_fontsize=8,
# )
# ax2.set_ylim(0.87, 4)
ax1.set_ylim(0.87, 40)
# ax_b.yaxis.set_major_locator(ticks_0p5)
DS_OVERLAPPED_XLIM_K = (0.26, 45)
ax1.set_xlim(0.1, 45)
# ax2.set_xlim(0.26, 45)

# show_png_on_axis(ax2, 'presentations/window_viz.png', (0.75, 0.4), size=0.04, use_data_coords=False, color=COLOR_WINDOWED)

ax1.text(0.7, 0.7, '$\phi=0.02$', transform=ax1.transAxes)
# ax2.text(0.7, 0.6, '$\phi=0.11$', transform=ax2.transAxes)

common.add_exponential_index_indicator(ax1, -2, (0.3, 10), 'k', x_limits=(0.01, 0.4), )


common.save_fig(fig, f'{path}/periodic_effects.pdf', hide_metadata=True, dpi=300)

########## drift vs theta #############
import particle_linking.drift_vs_theta
fig, ax = plt.subplots(figsize=(5, 4))
# files = common.files_from_filenames('particle_linking/data', 'trajs_', 'sim_shear0.0_T296_theta*_nblobs162_dt0.01_tmax1000')
# print(files)
particle_linking.drift_vs_theta.go(
    ['sim_shear0.0_T296_theta0_RHS_nblobs162_dt0.01_tmax1000',
     'sim_shear0.0_T296_theta2_RHS_nblobs162_dt0.01_tmax1000',
     'sim_shear0.0_T296_theta4_RHS_nblobs162_dt0.01_tmax1000',
     'sim_shear0.0_T296_theta6_RHS_nblobs162_dt0.01_tmax1000',
     'sim_shear0.0_T296_theta8_RHS_nblobs162_dt0.01_tmax1000',
     'sim_shear0.0_T296_theta10_RHS_nblobs162_dt0.01_tmax1000',],
    ax,
)
common.save_fig(fig, f'{path}/drift_vs_theta.pdf', hide_metadata=True)



########## drift vs phi ##############
import particle_linking.drift_vs_phi
fig, ax = plt.subplots(figsize=(5, 4))
# files = common.files_from_filenames('particle_linking/data', 'trajs_', 'sim_hydro_*_L640_theta10')
# print(files)
files = ['sim_hydro_002_L640_theta10', 'sim_hydro_004_L640_theta10', 'sim_hydro_006_L640_theta10', 'sim_hydro_008_L640_theta10', 'sim_hydro_010_L640_theta10']
particle_linking.drift_vs_phi.go([], ax, eleanor=True)
ax.set_ylim(0.47, 1.21)
ax.set_xlim(0.0, 0.6)
common.save_fig(fig, f'{path}/drift_vs_phi.pdf', hide_metadata=True)
particle_linking.drift_vs_phi.go(files, ax)
# ax.set_ylim(0.47, 1.21)
common.save_fig(fig, f'{path}/drift_vs_phi2.pdf', hide_metadata=True)




########### two walls ##########

WALL_LONG_COLOR = 'royalblue'
WALL_SHORT_COLOR = 'cornflowerblue'
OPEN_LONG_COLOR = 'indianred'
OPEN_SHORT_COLOR = 'lightcoral'
THEORY_COLOR = 'indigo'
WALLTRAP_LONG_COLOR = 'seagreen'
WALLTRAP_SHORT_COLOR = 'mediumseagreen'

COLOR_TWOWALLS_SHORT = 'mediumturquoise'
COLOR_TWOWALLS_LONG = 'lightseagreen'
fig, ax = plt.subplots(figsize=(4, 3))
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
        # dict(
        #     file='ld_hydro_dpstokes_0.114_twowalls_L1280_t8h_32s',
        #     source='f_t64',
        #     msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
        #     plot_index=np.index_exp[:58],
        #     # color=COLOR_TWOWALLS_LONG,
        #     label='two walls, $t=64\mathrm{{s}}$'
        # ),
        # dict(
        #     file='ld_hydro_dpstokes_0.114_twowalls_L1280_t8h_32s',
        #     source='f_t256',
        #     msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
        #     plot_index=np.index_exp[:58],
        #     # color=COLOR_TWOWALLS_LONG,
        #     label='two walls, $t=256\mathrm{{s}}$'
        # ),
        # dict(
        #     file='ld_hydro_dpstokes_0.114_twowalls_L1280_t8h_32s',
        #     source='f_t1024',
        #     msd_file='ld_hydro_dpstokes_0.114_twowalls_L1280_t1h_1_unwrap_2d',
        #     plot_index=np.index_exp[:58],
        #     # color=COLOR_TWOWALLS_LONG,
        #     label='two walls, $t=1024\mathrm{{s}}$'
        # ),
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

presentations.presentations.show_png_on_axis(ax, 'presentations/side_view-two_walls-1p5sigma.drawio.png', (0.75, 0.46), size=0.05, data_coords=False, color=COLOR_TWOWALLS_LONG)

ax.text(0.7, 0.63, '$\phi=0.11$', transform=ax.transAxes)


common.save_fig(fig, f'{path}/twowalls.pdf', hide_metadata=True, dpi=300)


####### drift vs shear ########
fig, ax = plt.subplots(figsize=(5, 4))
import particle_linking.drift_vs_shear
files = common.files_from_filenames('particle_linking/data', 'trajs_', 'sim_shear*_T0_RHS_nblobs642')
particle_linking.drift_vs_shear.go(files, ax)
common.save_fig(fig, f'{path}/drift_vs_shear.pdf', hide_metadata=True)



####### trap and wall ######
fig1, ax1 = plt.subplots(1, 1, figsize=(3, 3))
fig2, ax2 = plt.subplots(1, 1, figsize=(3, 3))
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
presentations.presentations.show_png_on_axis(ax1, '/home/acarter/presentations/donev/figures/side_view-potential_and_wall-annotated.drawio.png',        (0.8, 0.6), size=0.05, color=WALLTRAP_LONG_COLOR, data_coords=False)
presentations.presentations.show_png_on_axis(ax2, '/home/acarter/presentations/donev/figures/side_view-potential_and_wall-higher-annotated.drawio.png', (0.8, 0.6), size=0.05, color=WALLTRAP_LONG_COLOR, data_coords=False)
# show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (0.88, 0.6), 0.04, color=OPEN_LONG_COLOR, use_data_coords=False)

common.save_fig(fig1, f'{path}/potentialandwall1.pdf', hide_metadata=True, dpi=300)
common.save_fig(fig2, f'{path}/potentialandwall2.pdf', hide_metadata=True, dpi=300)


############ collective diffusion above one wall ####################

fig, ax = plt.subplots(figsize=(5, 4))

ax.set_ylim(0.8, 20)
ax.set_xlim(2.1e-2, 3e1)
ax.semilogy()

visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='tab:green'
        ),
    ),
    ax = ax,
    plot_against_k=True,
    legend_fontsize=7,
    discrete_colors=True,
    show_twin_k_axis=False,
    allow_rescale_x=True
)

common.save_fig(fig, f'{path}/collective_diffusion_1.pdf', hide_metadata=True)
common.save_fig(fig, f'{path}/collective_diffusion_1.png', hide_metadata=True, dpi=300)


visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_open_0.114_L2560_t1h_1',
            source='f_t1',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[19:],
            color=OPEN_SHORT_COLOR,
            label='potential confinement, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_open_0.114_L2560',
            source='f_t1024',
            msd_file='ld_hydro_nbody_open_0.114_L2560_t1h_1_2d',
            plot_index=np.index_exp[:19],
            color=OPEN_LONG_COLOR,
            label='potential confinement, $t=1024\mathrm{{s}}$'
        ),
    ),
    ax = ax,
    plot_against_k=True,
    legend_fontsize=7,
    discrete_colors=True,
    show_twin_k_axis=False,
    allow_rescale_x=True
)
common.add_exponential_index_indicator(ax, -1, (1e-1, 1e1), 'k', x_limits=(1e-2, 2e-1), )
presentations.presentations.show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (1e1, 5), 0.04, color=OPEN_LONG_COLOR)


common.save_fig(fig, f'{path}/collective_diffusion_2.pdf', hide_metadata=True)
common.save_fig(fig, f'{path}/collective_diffusion_2.png', hide_metadata=True, dpi=300)



visualisation.Ds_overlapped_mult.go(
    (
        dict(
            file='ld_hydro_nbody_0.114_L2560_t1h_1', # 3600 running
            source='f_t1',
            plot_index=np.index_exp[23:],
            color=WALL_SHORT_COLOR,
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d', # 3600 running, use 2d
            label='single wall, $t=1\mathrm{{s}}$'
        ),
        dict(
            file='ld_hydro_nbody_0.114_L2560',
            source='f_t1024',
            msd_file='ld_hydro_nbody_0.114_L2560_t1h_1_2d', # 3600 running, use 2d
            color=WALL_LONG_COLOR,
            label='single wall, $t=1024\mathrm{{s}}$'
        ),
    ),
    ax = ax,
    plot_against_k=True,
    legend_fontsize=7,
    discrete_colors=True,
    show_twin_k_axis=False,
    allow_rescale_x=True
)
presentations.presentations.show_png_on_axis(ax, 'presentations/libmobility/side_view_2.png', (1e1, 2.5), 0.04, color=WALL_LONG_COLOR)
# show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (1e1, 2.5), 0.04, color=WALL_LONG_COLOR)


common.save_fig(fig, f'{path}/collective_diffusion_3.pdf', hide_metadata=True)
common.save_fig(fig, f'{path}/collective_diffusion_3.png', hide_metadata=True, dpi=300)