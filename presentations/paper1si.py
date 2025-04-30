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
from presentations.paper1 import DS_OVERLAPPED_XLIM, DS_OVERLAPPED_XLIM_K, DS_OVERLAPPED_XLIM_L_FOR_FKT, DS_OVERLAPPED_XLIM1280, DS_OVERLAPPED_XLIM_L_FOR_FKT1280, \
    DS_OVERLAPPED_YLIM, DS_OVERLAPPED_YLIM_FKT, DS_OVERLAPPED_YLIM_FKT_SMALL, \
    COLOR_PHI011, COLOR_PHI011_THEORY, COLOR_PHI002, COLOR_PHI011_alt, \
    MARKER_COUNTING, MARKER_PHENOMFIT, MARKER_FKT, \
    LINESTYLE_COUNTING, \
    set_title_sim, set_title_exp, LEGEND_FONTSIZE

COLOR_LX_320 = cmocean.cm.ice(0.75)
COLOR_LX_640 = cmocean.cm.ice(0.5)
COLOR_LX_1280 = cmocean.cm.ice(0.25)

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
# file_labels = ['{0:.2g} hours'.format(common.load(f'visualisation/data/Ds_from_f_first_first_{file}.npz')['max_time_hours']) for file in files]

# visualisation.Ds_overlapped_mult.go(
#     files,
#     sources=['f_first_first'],
#     ax=ax,
#     logarithmic_y=True,
#     plot_against_k=True,
#     legend_fontsize=7,
#     # linestyle='-',
#     file_labels=file_labels,
#     errorbar_alpha=0,
#     discrete_colors=False,
# )
# common.save_fig(fig, f'{path}/si_f_divergence_movie_length.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_f_divergence_movie_length.pdf', hide_metadata=True)

# # ##########################f(k, t) short vs long #########################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# # visualisation.Ds_overlapped_mult.go(
# #     ['sim_nohydro_011_L320', 'sim_nohydro_011_L320_longer'],
# #     ax,
# #     ['f_first_first', 'f_long'],
# #     # PLOT_AGAINST_K=True,
# #     logarithmic_y=False,
# #     plot_against_k=True,
# # )
# #visualisation/data/Ds_from_f_long_sim_nohydro_011_L640_longer_mergedD.npz
# #visualisation/data/Ds_from_f_long_sim_nohydro_011_L320_longer_mergedD.npz
# visualisation.Ds_overlapped_mult.go(
#     # files = ['sim_nohydro_011_L640', 'sim_nohydro_011_L640_longer'],
#     files = ['sim_nohydro_011_L640_longer_mergedD'],
#     sources = ['f_first_first', 'f_long'],
#     ax=ax,
#     # PLOT_AGAINST_K=True,
#     logarithmic_y=False,
#     plot_against_k=True,
#     colors=[[COLOR_PHI011_alt, 'tab:orange']],#, ['tab:blue', 'tab:green']],
#     markers=MARKER_FKT,
#     source_labels=['$f(k, t)$ short time', '$f(k, t)$ long time'],
#     file_labels = [''],#, 'longer'],
# )
# ax.set_ylim(*DS_OVERLAPPED_YLIM_FKT_SMALL)
# common.save_fig(fig, f'{path}/si_fkt_short_long.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_fkt_short_long.pdf', hide_metadata=True)

############################################################
# # # # fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# # # # box_counting.D_of_L.go(
# # # #     'eleanorlong010',
# # # #     # 'var',
# # # #     'nmsdfitinter',
# # # #     ax=ax,
# # # #     save_data=False,
# # # #     title='',
# # # #     show_nofit_cutoff=False,
# # # #     max_num_boxes=6,
# # # # )

# # # # common.save_fig(fig, f'{path}/si_tsi_fits.png', hide_metadata=True)
# # # # common.save_fig(fig, f'{path}/si_tsi_fits.pdf', hide_metadata=True)


########################## high density periodic #################################
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
# visualisation.Ds_overlapped_mult.go(
#     # ['sim_nohydro_011_L320', 'sim_nohydro_011_L640_crop320'],
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
# fig, axs = plt.subplots(2, len(sources), figsize=(len(sources)*2.5, 2*3), sharey='row')
# box_counting.plateau_sources.go(
#     # file                 = 'eleanorlong010',
#     file                 = 'sim_nohydro_011_L1280', # first is merged then mergedD
#     axs                  = axs,
#     sources              = sources,
#     Ds_overlapped_kwargs = dict(
#         colors        = [ # sources are tsi_replace, tsi, theory
#                             # [COLOR_PHI011_alt, 'none',       'none'],
#                             # ['none',           COLOR_PHI011, COLOR_PHI011_THEORY],
#                             [COLOR_PHI011_alt, 'none',       'none'],
#                             ['none',           COLOR_PHI011, COLOR_PHI011_THEORY],
#                         ],
#         # file_labels   = [''],
#         source_labels = ['phenomelogical fit', 'timescale integral', 'theory'],
#         markers       = [[MARKER_PHENOMFIT, MARKER_COUNTING, 'none']]*2,
#         linestyles    = [['none', 'none', LINESTYLE_COUNTING]]*2,
#     ),
#     D_of_L_kwargs = dict(
#         late_C_N_alpha = 1,
#         show_long_fits = False,
#     ),
# )
# for ax in axs[1]:
#     ax.set_ylim(DS_OVERLAPPED_YLIM)
#     ax.set_xlim(DS_OVERLAPPED_XLIM)
# common.save_fig(fig, f'{path}/si_plateau_sources2.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateau_sources2.pdf', hide_metadata=True)


############### plateau vs L ###############
# """
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
# """


############## plateaus rescaled ####################
fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
for Lx, color in [
    [1280, COLOR_LX_1280],
    [640,  COLOR_LX_640],
    [320,  COLOR_LX_320],
]:
    box_counting.plateaus.go(
        f'sim_nohydro_011_L{Lx}_longer_frac_of_window',
        ax,
        ['var', 'sDFT'] if Lx==1280 else ['var'],
        rescale_window=True,
        label_prefix=f'$L_x = {Lx}\mathrm{{\mu m}}$ ',
        rescale_sigma=False,
        colors=[color, 'black'] if Lx==1280 else [color],
    )
ax.legend(fontsize=LEGEND_FONTSIZE)
ax.set_xlim(2e-2, 1e0)
ax.set_ylim(1e-5, 2e-2)
ax.set_title('periodic nonhydro')
common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled.pdf', hide_metadata=True)

# nonperiodic version
fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
for Lx, color in [
    # [1280, COLOR_LX_1280],
    [640,  COLOR_LX_640],
    [320,  COLOR_LX_320],
]:
    box_counting.plateaus.go(
        f'sim_nohydro_011_L{Lx}_pot_longer_frac_of_window',
        ax,
        ['var', 'sDFT'] if Lx==640 else ['var'],
        rescale_window=True,
        label_prefix=f'$L_x = {Lx}\mathrm{{\mu m}}$ ',
        rescale_sigma=False,
        colors=[color, 'black'] if Lx==640 else [color],
    )
ax.legend(fontsize=LEGEND_FONTSIZE)
ax.set_xlim(2e-2, 1e0)
ax.set_ylim(1e-5, 2e-2)
ax.set_title('nonperiodic nonhydro')
common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled_nonperiodic.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled_nonperiodic.pdf', hide_metadata=True)

# hydro periodic version
fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
for Lx, color in [
    # [1280, COLOR_LX_1280],
    [640,  COLOR_LX_640],
    [320,  COLOR_LX_320],
]:
    box_counting.plateaus.go(
        f'sim_hydro_011_L{Lx}_longer_frac_of_window',
        ax,
        ['var', 'sDFT'] if Lx==640 else ['var'],
        rescale_window=True,
        label_prefix=f'$L_x = {Lx}\mathrm{{\mu m}}$ ',
        rescale_sigma=False,
        colors=[color, 'black'] if Lx==640 else [color],
    )
ax.legend(fontsize=LEGEND_FONTSIZE)
ax.set_xlim(2e-2, 1e0)
ax.set_ylim(1e-5, 2e-2)
ax.set_title('periodic hydro')
common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled_periodic_hydro.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled_periodic.pdf', hide_metadata=True)

# nonperiodic hydro version
fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)
for Lx, color in [
    # [1280, COLOR_LX_1280],
    [640,  COLOR_LX_640],
    [320,  COLOR_LX_320],
]:
    box_counting.plateaus.go(
        f'sim_hydro_011_L{Lx}_pot_longer_frac_of_window',
        ax,
        ['var', 'sDFT'] if Lx==640 else ['var'],
        rescale_window=True,
        label_prefix=f'$L_x = {Lx}\mathrm{{\mu m}}$ ',
        rescale_sigma=False,
        colors=[color, 'black'] if Lx==640 else [color],
    )
ax.legend(fontsize=LEGEND_FONTSIZE)
ax.set_xlim(2e-2, 1e0)
ax.set_ylim(1e-5, 2e-2)
ax.set_title('nonperiodic hydro')
common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled_nonperiodic_hydro.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_plateaus_sim_rescaled_nonperiodic.pdf', hide_metadata=True)

ADD crop-.9_frac_of_window

########################## periodic box size ##################
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=SI_DOUBLE_FIGSIZE, sharey='row')

# # files = ['sim_nohydro_011_L160_longer_mergedD', 'sim_nohydro_011_L320_longer_mergedD', 'sim_nohydro_011_L640_longer_mergedD', 'sim_nohydro_011_L1280_longer_mergedD'],
# # Ls = [160, 320, 640, 1280]
# Ls = [320, 640, 1280]

# kwargs = dict(
#     file_labels = [f'$L_x={L}\mathrm{{\mu m}}$' for L in Ls],
#     # colors = [[cmocean.cm.ice((len(Ls)-i)/len(Ls))] for i in range(len(Ls))]
#     colors = [[COLOR_LX_320], [COLOR_LX_640], [COLOR_LX_1280]]#, [cmocean.cm.ice(0.1)]],
# )
# xlim = (10, 100)
# ylim = (2e1, 2e3)

# visualisation.Ds_overlapped_mult.go(
#     files   = [f'sim_nohydro_011_L{L}_longer_mergedD' for L in Ls],
#     sources = ['timescaleint_nofit_cropped_var'],
#     ax      = ax1,
#     markers = MARKER_COUNTING,
#     **kwargs
# )
# ax1.set_ylim(DS_OVERLAPPED_YLIM)
# ax1.set_xlim(DS_OVERLAPPED_XLIM1280)
# ax_label(ax1, 'a', x=0.95)
# # ax_sim.set_xlim(plateaus_xlim)
# # ax_sim.set_ylim(plateaus_ylim)
# # ax_sim.tick_params(axis='x', which='minor', labelbottom=False)

# visualisation.Ds_overlapped_mult.go(
#     # files   = [f'sim_nohydro_011_L{L}_mixt' for L in Ls],
#     files   = [f'sim_nohydro_011_L{L}_longer_mergedD' for L in Ls],
#     sources = ['f_first_first'],
#     ax      = ax2,
#     markers = MARKER_FKT,
#     disable_ylabel=True,
#     **kwargs
# )
# ax2.set_xlim(DS_OVERLAPPED_XLIM_L_FOR_FKT1280)
# ax2.set_ylim(DS_OVERLAPPED_YLIM_FKT)
# ax_label(ax2, 'b', x=0.95)
# # ax_exp.set_xlim(plateaus_xlim)
# # ax_exp.set_ylim(plateaus_ylim)
# # ax_exp.tick_params(axis='x', which='minor', labelbottom=False)

# common.save_fig(fig, f'{path}/si_D_periodic_size.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_D_periodic_size.pdf', hide_metadata=True)

##### nonperiodic version
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=SI_DOUBLE_FIGSIZE, sharey='row')

Ls = [640,]

kwargs = dict(
    labels = [f'$L_x={L}\mathrm{{\mu m}}$' for L in Ls],
    # colors = [[cmocean.cm.ice((len(Ls)-i)/len(Ls))] for i in range(len(Ls))]
    # colors = [[COLOR_LX_320], [COLOR_LX_640], [COLOR_LX_1280]]#, [cmocean.cm.ice(0.1)]],
    colors = [COLOR_LX_640]#, [cmocean.cm.ice(0.1)]],
)
xlim = (10, 100)
ylim = (2e1, 2e3)

visualisation.Ds_overlapped_mult.go(
    [(f'sim_nohydro_011_L{L}_pot_longer_mergedD', 'timescaleint_nofit_cropped_var') for L in Ls],
    ax      = ax1,
    markers = MARKER_COUNTING,
    **kwargs
)
ax1.set_ylim(DS_OVERLAPPED_YLIM)
ax1.set_xlim(DS_OVERLAPPED_XLIM1280)
ax_label(ax1, 'a', x=0.95)
# ax_sim.set_xlim(plateaus_xlim)
# ax_sim.set_ylim(plateaus_ylim)
# ax_sim.tick_params(axis='x', which='minor', labelbottom=False)

# visualisation.Ds_overlapped_mult.go(
#     # files   = [f'sim_nohydro_011_L{L}_mixt' for L in Ls],
#     files   = [f'sim_nohydro_011_L{L}_pot_longer_mergedD' for L in Ls],
#     sources = ['f_first_first'],
#     [(f'sim_nohydro_011_L{L}_pot_longer_mergedD', 'f_first_first') for L in Ls]
#     ax      = ax2,
#     markers = MARKER_FKT,
#     disable_ylabel=True,
#     **kwargs
# )
ax2.set_xlim(DS_OVERLAPPED_XLIM_L_FOR_FKT1280)
ax2.set_ylim(DS_OVERLAPPED_YLIM_FKT)
ax_label(ax2, 'b', x=0.95)
# ax_exp.set_xlim(plateaus_xlim)
# ax_exp.set_ylim(plateaus_ylim)
# ax_exp.tick_params(axis='x', which='minor', labelbottom=False)

common.save_fig(fig, f'{path}/si_D_periodic_size_nonperiodic.png', hide_metadata=True)
common.save_fig(fig, f'{path}/si_D_periodic_size_nonperiodic.pdf', hide_metadata=True)


# ########################## timescaleint longtime fit ##################
# # # # fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# # # # visualisation.Ds_overlapped_mult.go(
# # # #     files = ['eleanorlong001', 'eleanorlong010'],
# # # #     sources = ['timescaleint_nofit_cropped_var', 'timescaleint_fixexponent_var'],
# # # #     source_labels = ['standard integration', 'integration with long-time extension'],
# # # #     file_labels=['$\phi=0.02$', '$\phi=0.11$'],
# # # #     colors = [[COLOR_PHI002, COLOR_PHI002], [COLOR_PHI011, COLOR_PHI011]],
# # # #     ax      = ax,
# # # #     markers = [[MARKER_COUNTING, 'x']]*2,
# # # #     legend_fontsize=7
# # # # )
# # # # ax.set_ylim(DS_OVERLAPPED_YLIM)
# # # # ax.set_xlim(DS_OVERLAPPED_XLIM)


# # # # common.save_fig(fig, f'{path}/si_timescaleint_extension.png', hide_metadata=True)
# # # # common.save_fig(fig, f'{path}/si_timescaleint_extension.pdf', hide_metadata=True)

# ###################################################
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=SI_DOUBLE_FIGSIZE)#, sharey='row')

# args = dict(
#     rescale_x=box_counting.msd_single.RESCALE_X_L2,
#     rescale_y=box_counting.msd_single.RESCALE_Y_L2
# )
# # ylim = (1e-3*0.3, 0.7e-1*0.3)
# # xlim = (1e-3, 1e5)
# text_pos = (0.6, 0.1)

# box_counting.msd_single.go(
#     'eleanorlong001',
#     ax1,
#     **args
# )
# ax1.set_ylim(1e-3*0.3, 0.7e-1*0.3)
# ax1.set_xlim(1e-3, 1e5)
# ax_label(ax1, 'a', x=0.95)
# ax1.text(*text_pos,'$\phi=0.02$', fontsize=12, transform=ax1.transAxes)

# box_counting.msd_single.go(
#     'eleanorlong010',
#     ax2,
#     **args
# )
# ax2.set_ylim(2e-3, 1e-1)
# ax2.set_xlim(1e-3, 1e5)
# ax_label(ax2, 'b', x=0.95)
# ax2.text(*text_pos, '$\phi=0.11$', fontsize=12, transform=ax2.transAxes)
# # ax2.set_y4label(None)

# common.save_fig(fig, f'{path}/si_nmsd_rescaled.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_nmsd_rescaled.pdf', hide_metadata=True)


########### hydro ############
# fig, ax = plt.subplots(1, 1, figsize=SI_SINGLE_FIGSIZE)

# ### need to do show_f_first here

# visualisation.Ds_overlapped_mult.go(
#     files = ['sim_nohydro_011_L640_longer_mergedD', 'sim_hydro_011_L640'],
#     sources = ['f_first_first'],
#     # source_labels = ['standard integration', 'integration with long-time extension'],
#     file_labels=['no hydro', 'hydro'],
#     # colors = [[COLOR_PHI002, COLOR_PHI002], [COLOR_PHI011, COLOR_PHI011]],
#     ax      = ax,
#     # markers = [[MARKER_COUNTING, 'x']]*2,
#     legend_fontsize=10
# )
# ax.set_ylim(DS_OVERLAPPED_YLIM)
# ax.set_xlim(DS_OVERLAPPED_XLIM_L_FOR_FKT)

# common.save_fig(fig, f'{path}/si_hydro.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/si_hydro.pdf', hide_metadata=True)