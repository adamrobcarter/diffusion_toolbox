import common

import isf.show_Fs_overlayed
import visualisation.Dc_literature
import visualisation.Ds_overlapped
import visualisation.Ds_overlapped_mult
import box_counting.D_of_L
import box_counting.msd_single
import box_counting.msd_combined
import box_counting.quantify_overlap_show
import box_counting.show_raw_counts
import box_counting.plateau_sources
import box_counting.plateaus

import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec
import cmocean


path = 'presentations/paper1'
w = 12
plt.rcParams.update({'axes.labelsize': 12})

# DS_OVERLAPPED_YLIM = (0.8, 2.5)
# DS_OVERLAPPED_XLIM = (0.1, 50)
from presentations.paper1 import DS_OVERLAPPED_XLIM, DS_OVERLAPPED_XLIM_K, DS_OVERLAPPED_YLIM, DS_OVERLAPPED_YLIM_FKT, \
    COLOR_PHI011, COLOR_PHI011_THEORY, MARKER_COUNTING, MARKER_PHENOMFIT, \
    set_title_sim, set_title_exp, LEGEND_FONTSIZE

def ax_label(ax, label, x=0.05, y=0.95, color='black'):
    if x < 0.5:
        ha = 'left'
    elif x == 0.5:
        ha = 'center'
    else:
        ha = 'right'

    ax.text(
        x, y, label, transform=ax.transAxes,
        # add to transform + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans) to move in px
        fontsize=25, va='top', ha=ha, fontfamily='serif', zorder=100, color=color
    )

def show_png(ax, file):
    image = plt.imread(file)
    print(image.shape)
    ax.imshow(image, interpolation='none')
    ax.axis('off') 

SI_SINGLE_FIGSIZE = (4, 3.5)
SI_DOUBLE_FIGSIZE = (6.5, 3.5)

# ######################## f divergence movie length ###########################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
# files = ['eleanorlong001_trim1.0', 'eleanorlong001_trim0.5', 'eleanorlong001_trim0.25', 'eleanorlong001_trim0.125', 'eleanorlong001_trim0.0625']
# times = [str(common.load(f'visualisation/data/Ds_from_f_first_first_{file}.npz')['max_time_hours']) + ' hours' for file in files]

# visualisation.Ds_overlapped_mult.go(
#     files,
#     sources=['MSD_first', 'f_first_first'],
#     ax=ax,
#     logarithmic_y=True,
#     plot_against_k=True,
#     legend_fontsize=7,
#     # linestyle='-',
#     file_labels=times,
#     errorbar_alpha=0,
# )
# common.save_fig(fig, f'{path}/si_f_divergence_movie_length.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_f_divergence_movie_length.pdf', hide_metadata=True)

# # ###################################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# box_counting.msd_single.go(
#     'eleanorlong001',
#     ax,
#     rescale_x=box_counting.msd_single.RESCALE_X_L2,
#     rescale_y=box_counting.msd_single.RESCALE_Y_L2,
#     show_rescaled_theory=True,
#     # max_boxes_on_plot=6,
#     # box_size_indices=[0,1,2,3,4,5,6,7,8]
# )
# ax.set_ylim(1e-3*0.3, 0.7e-1*0.3)
# ax.set_xlim(1e-3, 1e5)
# common.save_fig(fig, f'{path}/si_nmsd_lowdensity.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_nmsd_lowdensity.pdf', hide_metadata=True)

# # ###################################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# box_counting.msd_single.go(
#     'eleanorlong010',
#     ax,
#     rescale_x=box_counting.msd_single.RESCALE_X_L2,
#     rescale_y=box_counting.msd_single.RESCALE_Y_L2
# )
# ax.set_ylim(1e-3, 0.7e-1)
# ax.set_xlim(1e-3, 1e5)
# common.save_fig(fig, f'{path}/si_nmsd_highdensity.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_nmsd_highdensity.pdf', hide_metadata=True)

# # ###################################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# # visualisation.Ds_overlapped_mult.go(
# #     ['sim_nohydro_011_L320', 'sim_nohydro_011_L320_longer'],
# #     ax,
# #     ['f_first_first', 'f_long'],
# #     # PLOT_AGAINST_K=True,
# #     logarithmic_y=False,
# #     plot_against_k=True,
# # )
# visualisation.Ds_overlapped.go(
#     'sim_nohydro_011_L320',
#     ['f_first_first', 'f_long'],
#     ax=ax,
#     # PLOT_AGAINST_K=True,
#     logarithmic_y=False,
#     plot_against_k=True,
# )
# common.save_fig(fig, f'{path}/si_fkt_short_long.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_fkt_short_long.pdf', hide_metadata=True)

# ############################################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# box_counting.D_of_L.go(
#     'eleanorlong010',
#     # 'var',
#     'nmsdfitinter',
#     ax=ax,
#     save_data=False,
#     title='',
#     show_nofit_cutoff=False,
#     max_num_boxes=6,
# )
# common.save_fig(fig, f'{path}/si_tsi_fits.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_tsi_fits.pdf', hide_metadata=True)

########################### no hydro f periodic #################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# visualisation.Ds_overlapped_mult.go(
#     ['sim_nohydro_002_L320', 'sim_nohydro_002_L640_crop320'],
#     ax,
#     ['MSD_first', 'f_first_first',],
#     logarithmic_y=False,
#     plot_against_k=True,
    # colors=[['olivedrab'], ['darkgreen']]
    # file_labels=['periodic boundary', 'non periodic boundary'],
# )
# ax.set_ylim(*DS_OVERLAPPED_YLIM_FKT)
# common.save_fig(fig, f'{path}/si_no_hydro_f_periodic.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_no_hydro_f_periodic.pdf', hide_metadata=True)


########################### high density periodic #################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
# visualisation.Ds_overlapped_mult.go(
#     ['sim_nohydro_011_L320', 'sim_nohydro_011_L640_crop320'],
#     sources=[
#         'f_first_first'
#     ],
#     ax=ax,
#     plot_against_k=True,
#     # markers=MARKER_FKT,
#     file_labels=['periodic boundary', 'non periodic boundary'],
#     source_labels=[''],
#     colors=[['olivedrab'], ['darkgreen']]
# )
# ax.set_ylim(*DS_OVERLAPPED_YLIM)
# # ax.yaxis.set_major_locator(ticks_0p5)
# ax.set_xlim(*DS_OVERLAPPED_XLIM_K)
# ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1, 0.8)) # puts upper right at specified position
# common.save_fig(fig, f'{path}/si_no_hydro_f_periodic_highdensity.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_no_hydro_f_periodic_highdensity.pdf', hide_metadata=True)


#################### plateau sources #####################
# sources = ['var', 'varmod', 'nmsdfit', 'sDFT']
# fig, axs = plt.subplots(2, len(sources), figsize=(len(sources)*3, 2*3))
# box_counting.plateau_sources.go(
#     file                 = 'eleanorlong010',
#     axs                  = axs,
#     sources              = sources,
#     Ds_overlapped_kwargs = dict(
#         colors        = [[COLOR_PHI011, 'black', COLOR_PHI011, COLOR_PHI011_THEORY]],
#         file_labels   = [''],
#         source_labels = ['phenom. fit', 'tsi nofit', 'tsi fix exp', 'theory'],
#         markers       = [[MARKER_PHENOMFIT, MARKER_COUNTING, MARKER_COUNTING, 'none']]
#     ),
# )
# for ax in axs[1]:
#     ax.set_ylim(DS_OVERLAPPED_YLIM)
#     ax.set_xlim(DS_OVERLAPPED_XLIM)
# common.save_fig(fig, f'{path}/si_plateau_sources.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateau_sources.pdf', hide_metadata=True)


############### plateau vs L ###############
# fig, (ax_sim, ax_exp) = plt.subplots(1, 2, figsize=SI_DOUBLE_FIGSIZE)

# plateaus_kwargs = dict(
#     sources = ['var', 'sDFT']
# )
# plateaus_xlim = (10, 100)
# plateaus_ylim = (2e1, 2e3)

# box_counting.plateaus.go(
#     'sim_nohydro_011_L640_longer',
#     ax_sim,
#     **plateaus_kwargs
# )
# set_title_sim(ax_sim)
# ax_sim.set_xlim(plateaus_xlim)
# ax_sim.set_ylim(plateaus_ylim)
# ax_sim.tick_params(axis='x', which='minor', labelbottom=False)

# box_counting.plateaus.go(
#     'eleanorlong010',
#     ax_exp,
#     **plateaus_kwargs
# )
# set_title_exp(ax_exp)
# ax_exp.set_xlim(plateaus_xlim)
# ax_exp.set_ylim(plateaus_ylim)
# ax_exp.tick_params(axis='x', which='minor', labelbottom=False)

# common.save_fig(fig, f'{path}/si_plateaus.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateaus.pdf', hide_metadata=True)


############### plateaus rescaled ####################
fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
box_counting.plateaus.go(
    'sim_nohydro_011_L640_longer_frac_of_window',
    ax,
    ['var', 'sDFT'],
    rescale_window=True,
    label_prefix='$L_x = 640\mathrm{\mu m}$ ',
    rescale_sigma=False,
    colors=[cmocean.cm.ice(0.25), cmocean.cm.ice(0.25)]
)
box_counting.plateaus.go(
    'sim_nohydro_011_L320_longer_frac_of_window',
    ax,
    ['var', 'sDFT'],
    rescale_window=True,
    label_prefix='$L_x = 320\mathrm{\mu m}$ ',
    rescale_sigma=False,
    colors=[cmocean.cm.ice(0.75), cmocean.cm.ice(0.75)],
    linestyle='none'
)
ax.legend(fontsize=LEGEND_FONTSIZE)
ax.set_xlim(2e-2, 1e0)
ax.set_ylim(1e-5, 2e-2)
common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled.png', hide_metadata=True)
common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled.pdf', hide_metadata=True)