import common
import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
# from presentations.paper1 import DS_OVERLAPPED_YLIM, DS_OVERLAPPED_XLIM, DS_OVERLAPPED_XLIM_K, DS_OVERLAPPED_YLIM_FKT, \
#     MARKER_COUNTING, MARKER_FKT, COLOR_PHI002, COLOR_PHI002_THEORY, COLOR_PHI011, COLOR_PHI011_THEORY, COLOR_PHI011_alt
from presentations.paper1 import *
import particle_detection.show

path = '/home/acarter/presentations/cmif25/figures'

WIDTH *= 1.2
ONEROW_HEIGHT *= 1.15

COLOR_PHI011_COUNTING        = common.colormap(0.3)
COLOR_PHI011_COUNTING_THEORY = 'black'
COLOR_PHI002_COUNTING        = common.colormap(0.7)
COLOR_PHI002_COUNTING_THEORY = 'black'
COLOR_PHI011_FKT             = common.colormap_cool(0.3)
COLOR_PHI011_FKT_THEORY      = 'black'
COLOR_PHI002_FKT             = common.colormap_cool(0.8)
COLOR_PHI002_FKT_THEORY      = 'black'


# fig, ax = plt.subplots(1, 1, figsize=(5, 4))
# # Ds_ov_mult_props = dict(
# #     logarithmic_y=False,
# #     legend_fontsize=LEGEND_FONTSIZE,
# #     file_labels=PHI_FILE_LABELS,
# #     source_labels=['theory', 'eq. 8'],
# #     colors=[[COLOR_PHI002_THEORY, COLOR_PHI002], [COLOR_PHI011_THEORY, COLOR_PHI011]],
# #     markers=[['none', MARKER_COUNTING]] * 2,
# #     linestyles=[[LINESTYLE_COUNTING, 'none']]*2,
# #     sources=['D_of_L_theory', 'timescaleint_nofit_cropped_var'],
# # )
# visualisation.Ds_overlapped_mult.go(
#     # ['eleanorlong001', 'eleanorlong010'],
#     [
#         ('eleanorlong010', 'D_of_L_theory'),
#         ('eleanorlong010', 'boxcounting_first_quad'),
#         ('eleanorlong010', 'timescaleint_nofit_cropped_var'),
#         ('sim_nohydro_011_L1280_longer_mergedD', 'timescaleint_nofit_cropped_var')
#     ],
#     ax=ax,
#     labels=['theory', 'experiment', 'experiment', 'simulation (no hydro)'],
#     discrete_colors=True,
#     # **Ds_ov_mult_props,
# )
# ax.set_ylim(*DS_OVERLAPPED_YLIM)
# # ax_b.yaxis.set_major_locator(ticks_0p5)
# ax.set_xlim(*DS_OVERLAPPED_XLIM)
# # ax_b.xaxis.labelpad = -1
# # ax_label(ax_b, 'd', x=0.95)
# # # ax_label(ax_b, 'd')
# # ax_b.get_legend().set_bbox_to_anchor((-0.05, 0, 1, 1))
# # ax_b.get_legend().set_bbox_to_anchor((-0.9, 0.01, 1, 1))
# common.save_fig(fig, f'{path}/tsi.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/tsi.pdf', hide_metadata=True)

# ################## fkt plot ##################
# fig, ax = plt.subplots(1, 1, figsize=(5, 4))
# # Ds_ov_mult_props = dict(
# #     logarithmic_y=False,
# #     legend_fontsize=LEGEND_FONTSIZE,
# #     file_labels=PHI_FILE_LABELS,
# #     source_labels=['theory', 'eq. 8'],
# #     colors=[[COLOR_PHI002_THEORY, COLOR_PHI002], [COLOR_PHI011_THEORY, COLOR_PHI011]],
# #     markers=[['none', MARKER_COUNTING]] * 2,
# #     linestyles=[[LINESTYLE_COUNTING, 'none']]*2,
# #     sources=['D_of_L_theory', 'timescaleint_nofit_cropped_var'],
# # )
# visualisation.Ds_overlapped_mult.go(
#     # ['eleanorlong001', 'eleanorlong010'],
#     [
#         ('eleanorlong010', 'D0Sk_theory'),
#         ('eleanorlong010', 'f_first_first'),
#         ('sim_nohydro_011_L1280_longer_mergedD', 'f_first_first')
#     ],
#     ax=ax,
#     # colors=[COLOR_PHI011_THEORY, COLOR_PHI011_alt, COLOR_PHI011_alt]
#     # **Ds_ov_mult_props,
# )
# ax.set_ylim(*DS_OVERLAPPED_YLIM_FKT)
# # ax_b.yaxis.set_major_locator(ticks_0p5)
# ax.set_xlim(*DS_OVERLAPPED_XLIM_K)
# # ax_b.xaxis.labelpad = -1
# # ax_b.get_legend().set_bbox_to_anchor((-0.05, 0, 1, 1))
# # ax_b.get_legend().set_bbox_to_anchor((-0.9, 0.01, 1, 1))
# common.save_fig(fig, f'{path}/fkt.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fkt.pdf', hide_metadata=True)



######## small windowing plot ##########
# fig, ax = plt.subplots(1, 1, figsize=(3, 3))
# # Ds_ov_mult_props = dict(
# #     logarithmic_y=False,
# #     legend_fontsize=LEGEND_FONTSIZE,
# #     file_labels=PHI_FILE_LABELS,
# #     source_labels=['theory', 'eq. 8'],
# #     colors=[[COLOR_PHI002_THEORY, COLOR_PHI002], [COLOR_PHI011_THEORY, COLOR_PHI011]],
# #     markers=[['none', MARKER_COUNTING]] * 2,
# #     linestyles=[[LINESTYLE_COUNTING, 'none']]*2,
# #     sources=['D_of_L_theory', 'timescaleint_nofit_cropped_var'],
# # )
# visualisation.Ds_overlapped_mult.go(
#     # ['eleanorlong001', 'eleanorlong010'],
#     [
#         ('eleanorlong010', 'D0Sk_theory'),
#         ('eleanorlong010', 'f_first_first'),
#         ('eleanorlong010_bhwindow', 'f_first_first'),
#         # ('eleanorlong010_mirrortile', 'f_first_first'),
#     ],
#     labels=['theory', 'exp.', 'exp., BH window'],
#     ax=ax,
#     # **Ds_ov_mult_props,
#     markers=['none', MARKER_FKT, MARKER_FKT],
#     colors=[COLOR_PHI011_FKT_THEORY, COLOR_PHI011_FKT, 'tab:blue'],
#     legend_fontsize=9,
# )
# ax.set_ylim(0.87, 4)
# # ax_b.yaxis.set_major_locator(ticks_0p5)
# ax.set_xlim(*DS_OVERLAPPED_XLIM_K)
# # ax_b.xaxis.labelpad = -1
# # ax_label(ax_b, 'd', x=0.95)
# # # ax_label(ax_b, 'd')
# # ax_b.get_legend().set_bbox_to_anchor((-0.05, 0, 1, 1))
# # ax_b.get_legend().set_bbox_to_anchor((-0.9, 0.01, 1, 1))
# common.save_fig(fig, f'{path}/windowing.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/windowing.pdf', hide_metadata=True)


########### comparison plot ###################

# fig, ax = plt.subplots(1, 1, figsize=(5, 4))
# visualisation.Ds_overlapped_mult.go(
#     # ['eleanorlong001', 'eleanorlong010'],
#     [
#         ('sim_nohydro_011_L1280_longer_mergedD', 'D_of_L_theory'),
#         ('eleanorlong010', 'timescaleint_nofit_cropped_var'),
#         ('sim_nohydro_011_L1280_longer_mergedD', 'timescaleint_nofit_cropped_var'),
#         ('sim_nohydro_011_L1280_longer_mergedD', 'D0Sk_theory'),
#         ('eleanorlong010', 'f_first_first'),
#         ('sim_nohydro_011_L1280_longer_mergedD', 'f_first_first'),
#     ],
#     ax=ax,
#     labels=['Countoscope theory', 'Countoscope experiment', 'Countoscope simulation', '$f(k, t)$ theory', '$f(k, t)$ experiment', '$f(k, t)$ simulation'],
#     colors=[COLOR_PHI011_THEORY, COLOR_PHI011, COLOR_PHI011_alt, COLOR_PHI011_THEORY, 'seagreen', 'limegreen'],
#     markers=[MARKER_COUNTING, MARKER_COUNTING, MARKER_COUNTING, MARKER_FKT, MARKER_FKT, MARKER_FKT],
#     # **Ds_ov_mult_props,
# )
# ax.set_ylim(*DS_OVERLAPPED_YLIM)
# ax.set_xlim(*DS_OVERLAPPED_XLIM)

# common.save_fig(fig, f'{path}/comparison.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/comparison.pdf', hide_metadata=True)

# ########################## fig tsi #########################
# fig, (ax_a, ax_b) = plt.subplots(1, 2,
#                     figsize=(WIDTH/2, ONEROW_HEIGHT),
#                     sharey='row',
#                     gridspec_kw={'wspace': GRIDSPEC_WSPACE, 'hspace': GRIDSPEC_HSPACE},
#                     )
# set_title_sim(ax_a)
# set_title_exp(ax_b)

# D_of_L_props = dict(
#     save_data=False,
#     title='',
#     show_nofit_cutoff=False,
#     max_num_boxes=6,
#     # colormap=COLORMAP_BOXES,
#     show_short_fits=True,
#     show_long_fits=False,
#     late_C_N_alpha=1,
# )

# # ax_a.set_xlim(*DS_OVERLAPPED_XLIM)
# phi = 0.114**3
# D0_over_Sk = (1 + phi) / ( 1- phi)
# # ax_a.hlines(D0_over_Sk, *ax_a.get_xlim(), color='gray', linestyle=':', label='$(1+\phi)/(1-\phi)^3$')

# # Ds_ov_mult_props = dict(
# #     logarithmic_y=False,
# #     legend_fontsize=LEGEND_FONTSIZE-1,
# #     file_labels=PHI_FILE_LABELS,
# #     source_labels=['theory', 'tsi no fit', 'tsi fix exponent'],
# #     colors=[[COLOR_PHI002_THEORY, COLOR_PHI002, 'black'], [COLOR_PHI011_THEORY, COLOR_PHI011, 'black']],
# #     markers=[['none', MARKER_COUNTING, 'x']] * 2,
# #     linestyles=[[LINESTYLE_COUNTING, 'none', 'none']]*2,
# #     sources=['D_of_L_theory', 'timescaleint_fixexponent_var', 'timescaleint_nofit_cropped_var'],
# # )
# Ds_ov_mult_props = dict(
#     logarithmic_y=False,
#     legend_fontsize=LEGEND_FONTSIZE,
#     labels=['$\phi=0.02$ theory', '$\phi=0.02$ observation', '$\phi=0.11$ theory', '$\phi=0.11$ observation'],
#     colors = [COLOR_PHI002_COUNTING_THEORY, COLOR_PHI002_COUNTING, COLOR_PHI011_COUNTING_THEORY, COLOR_PHI011_COUNTING],
#     markers=['none', MARKER_COUNTING]*2,
#     linestyles=[LINESTYLE_COUNTING, 'none']*2,
# )

# visualisation.Ds_overlapped_mult.go(
#     [
#         (f'{SIMULATION_SOURCE_002}_longer_{MERGE_TYPE_TIMESCALEINT}', 'D_of_L_theory'),
#         (f'{SIMULATION_SOURCE_002}_longer_{MERGE_TYPE_TIMESCALEINT}', 'timescaleint_nofit_cropped_var'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'D_of_L_theory'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'timescaleint_nofit_cropped_var'),
#     ],
#     ax=ax_a,
#     **Ds_ov_mult_props
#     # markers=[[None, MARKER_COUNTING, 'x'], [None, MARKER_COUNTING, 'x']]
# )
# ax_a.set_xlim(*DS_OVERLAPPED_XLIM)
# ax_a.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_a.yaxis.set_major_locator(ticks_0p5)
# # ax_b.yaxis.labelpad = -12
# ax_a.xaxis.labelpad = -1
# ax_a.get_legend().remove()

# visualisation.Ds_overlapped_mult.go(
#     [
#         ('eleanorlong001', 'D_of_L_theory'),
#         ('eleanorlong001', 'timescaleint_nofit_cropped_var'),
#         ('eleanorlong010', 'D_of_L_theory'),
#         ('eleanorlong010', 'timescaleint_nofit_cropped_var'),
#     ],
#     ax=ax_b,
#     **Ds_ov_mult_props,
#     disable_ylabel=True,
#     # markers=[[None, 'none', 'o'], [None, 'none', 'o'], [None, 'none', 'o']],
# )
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_b.yaxis.set_major_locator(ticks_0p5)
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
# ax_b.xaxis.labelpad = -1
# # ax_label(ax_b, 'd')
# ax_b.get_legend().set_bbox_to_anchor((-0.3, 0, 1, 1))
# # ax_b.get_legend().set_bbox_to_anchor((-0.9, 0.01, 1, 1))

# common.save_fig(fig, f'{path}/fig_tsi.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig_tsi.pdf', hide_metadata=True)




# ######################### fig 8 #########################
# fig, (ax_a, ax_b) = plt.subplots(1, 2,
#                                     figsize=(WIDTH/2, ONEROW_HEIGHT),
#                                     sharey='row',
#                                     gridspec_kw={'wspace': GRIDSPEC_WSPACE},
#                                     )
# set_title_sim(ax_a)
# set_title_exp(ax_b)

# fig8_props = dict(
#     logarithmic_y=False,
#     legend_fontsize=LEGEND_FONTSIZE,
#     # file_labels=[''],
#     # sources=      ['timescaleint_nmsdfitinter', 'D_of_L_theory',      'f_first_first', 'D0Sk_theory'],
#     # sources=      ['timescaleint_nofit_cropped_var',          'D_of_L_theory',      'f_first_first', 'D0Sk_theory'],
#     colors     = [COLOR_PHI011_COUNTING,     COLOR_PHI011_COUNTING_THEORY, COLOR_PHI011_FKT,        COLOR_PHI011_FKT_THEORY],
#     markers    = [MARKER_COUNTING,           'none',                       MARKER_FKT,              'none', ],
#     linestyles = ['none',                    LINESTYLE_COUNTING,           'none',                  LINESTYLE_FKT],
#     labels     = ['Countoscope observation', 'Countoscope theory',         '$f(k, t)$ observation', '$f(k, t)$ theory']
# )

# visualisation.Ds_overlapped_mult.go(
#     [
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'timescaleint_nofit_cropped_var'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'D_of_L_theory'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'f_first_first'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'D0Sk_theory'),
#     ],
#     ax=ax_a,
#     **fig8_props
# )
# ax_a.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_a.yaxis.set_major_locator(ticks_0p5)
# ax_a.set_xlim(*DS_OVERLAPPED_XLIM)
# ax_a.yaxis.labelpad = -1
# ax_a.xaxis.labelpad = -1
# ax_a.get_legend().remove()

# visualisation.Ds_overlapped_mult.go(
#     [
#         ('eleanorlong010', 'timescaleint_nofit_cropped_var'),
#         ('eleanorlong010', 'D_of_L_theory'),
#         ('eleanorlong010', 'f_first_first'),
#         ('eleanorlong010', 'D0Sk_theory'),
#     ],
#     ax=ax_b,
#     disable_ylabel=True,
#     **fig8_props
# )
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_b.yaxis.set_major_locator(ticks_0p5)
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
# ax_b.xaxis.labelpad = -1
# ax_b.get_legend().set_bbox_to_anchor((-0.4, 0, 1, 1))

# common.save_fig(fig, f'{path}/fig8.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig8.pdf', hide_metadata=True)



# # ######################## fig 7 (fkt) #########################
# fig, (ax_a, ax_b) = plt.subplots(1, 2,
#                                     figsize=(WIDTH/2, ONEROW_HEIGHT),
#                                     sharey='row',
#                                     gridspec_kw={'wspace': GRIDSPEC_WSPACE},
#                                     )
# set_title_sim(ax_a)
# set_title_exp(ax_b)

# fig7_props = dict(
#     logarithmic_y=False,
#     legend_fontsize=LEGEND_FONTSIZE,
#     colors=[COLOR_PHI002_FKT, COLOR_PHI002_FKT_THEORY, COLOR_PHI011_FKT, COLOR_PHI011_FKT_THEORY],
#     # file_labels=PHI_FILE_LABELS,
#     markers=[MARKER_FKT, 'none']*2,
#     linestyles=['none', LINESTYLE_FKT]*2,
#     labels=['$\phi=0.02$ observation', '$\phi=0.02$ theory', '$\phi=0.11$ observation', '$\phi=0.11$ theory'],
#     plot_against_k = True
# )

# visualisation.Ds_overlapped_mult.go(
#     [
#         (f'{SIMULATION_SOURCE_002}_longer_{MERGE_TYPE_TIMESCALEINT}', 'f_first_first'),
#         (f'{SIMULATION_SOURCE_002}_longer_{MERGE_TYPE_TIMESCALEINT}', 'D0Sk_theory'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'f_first_first'),
#         (f'{SIMULATION_SOURCE_011}_longer_{MERGE_TYPE_TIMESCALEINT}', 'D0Sk_theory'),
#     ],
#     ax=ax_a,
#     **fig7_props
# )
# ax_a.set_ylim(*DS_OVERLAPPED_YLIM_FKT_SMALL)
# ax_a.set_xlim(*DS_OVERLAPPED_XLIM_L_FOR_FKT)
# # ax_a.yaxis.labelpad = -5
# ax_a.xaxis.labelpad = -5
# ax_a.get_legend().remove()

# visualisation.Ds_overlapped_mult.go(
#     [
#         ('eleanorlong001', 'f_first_first'),
#         ('eleanorlong001', 'D0Sk_theory'),
#         ('eleanorlong010', 'f_first_first'),
#         ('eleanorlong010', 'D0Sk_theory'),
#     ],
#     ax=ax_b,
#     disable_ylabel=True,
#     **fig7_props
# )
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM_FKT_SMALL)
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM_K)
# # ax_b.yaxis.labelpad = -5
# ax_b.xaxis.labelpad = -5
# ax_b.get_legend().set_bbox_to_anchor((-0.85, 0, 1, 1))

# common.save_fig(fig, f'{path}/fig7.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig7.pdf', hide_metadata=True)


# ####### window viz #######
# fig, ax = plt.subplots(1, 1)

# particle_detection.show.go('eleanorlong010', ax, fig, bhwindow=True)

# common.save_fig(fig, f'{path}/window_viz.png', dpi=300, only_plot=True)


#### tracking video ###
# import particle_detection.show_movie
# particle_detection.show_movie.go(
#     'eleanor0.01',
#     infile = f'particle_linking/data/trajs_eleanor0.01.npz',
#     outfile = f'{path}/movie_eleanor0.01',
#     output_type='frames',
#     crop = 100,
#     every_nth_frame=10,
#     annotation_color='white',
# )
# particle_detection.show_movie.go(
#     'eleanor0.34',
#     infile = f'particle_linking/data/trajs_eleanor0.34.npz',
#     outfile = f'{path}/movie_eleanor0.34',
#     output_type='frames',
#     crop = 100,
#     every_nth_frame=10,
#     annotation_color='white',
# )


import preprocessing.stack_movie
data = common.load(f'preprocessing/data/stack_pierre_exp.npz')
preprocessing.stack_movie.go(
    'pierre_exp',
    data,
    outputfilename=f'{path}/movie_pierre_exp',
    output_type='frames',
)