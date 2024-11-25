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

import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec

"""
Figures should be prepared with the PDF layout in mind.
Individual figures should not be longer than one page and with a width that corresponds to
1 column (85 mm) or 2 columns (180 mm).

All images must have a resolution of 300 dpi at final size.
"""

path = 'presentations/paper1'
w = 12
plt.rcParams.update({'axes.labelsize': 12})

DS_OVERLAPPED_YLIM = (0.8, 2.5)
DS_OVERLAPPED_XLIM = (0.05, 50)
ticks_0p5 = plt.MultipleLocator(base=0.5)

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
    ax.imshow(image, interpolation='none')
    ax.axis('off') 
    # ax.axes.xaxis.set_visible(False) # you can use these two to hide the axes
    # ax.axes.yaxis.set_visible(False) # but keep a black border around the plot

def show_png2(fig, file):
    image = plt.imread(file)
    print(image.shape)
    # fig.figimage(image, xo=50, yo=50,  origin='lower')

    ax = fig.add_axes([0.02,0.05, 0.3, 0.45])#, frameon=True)
    ax.imshow(image, interpolation='none')
    ax_label(ax, 'c')
        # ax.axis('off') 
    # ax.imshow(image, interpolation='none')
    # ax.axis('off') 
    ax.axes.xaxis.set_visible(False) # you can use these two to hide the axes
    ax.axes.yaxis.set_visible(False) # but keep a black border around the plot


########################## fig 1 #########################
fig = plt.figure(figsize=(w/2, 5))
gs = matplotlib.gridspec.GridSpec(2, 1)#, height_ratios=(1, 1))#, wspace=0, hspace=0)
gs0 = gs[0].subgridspec(1, 3)#, width_ratios=[1, 1])
gs1 = gs[1].subgridspec(1, 2, width_ratios=[0.9, 1.1])
ax_a = fig.add_subplot(gs0[0])
ax_b = fig.add_subplot(gs0[1])
ax_c = fig.add_subplot(gs0[2])
ax_d = fig.add_subplot(gs1[0])
ax_e = fig.add_subplot(gs1[1])

show_png(ax_a, 'presentations/paper1/fig1_dhont1.png')
ax_label(ax_a, 'a', x=0.9)
show_png(ax_b, 'presentations/paper1/fig1_dhont2.png')
ax_label(ax_b, 'b', x=0.9)
show_png(ax_c, 'presentations/paper1/fig1_dhont3.png')
ax_label(ax_c, 'c', x=0.9)
show_png(ax_d, 'presentations/paper1/fig1_sketch.png')
ax_label(ax_d, 'd', x=0.9)
box_counting.show_raw_counts.go('eleanorlong010_trim0.0625', ax_e)
ax_label(ax_e, 'e')
ax_e.set_aspect(15)

# fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
# fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0, wspace=0)

common.save_fig(fig, f'{path}/fig1.png', hide_metadata=True)
common.save_fig(fig, f'{path}/fig1.pdf', hide_metadata=True)


######################### fig 2 #########################
# fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(w/2, 5.15),
#                                                  constrained_layout=True,
#                                                  gridspec_kw={'wspace': 0.09, 'hspace': 0})

# show_png(ax_a, 'presentations/paper1/image_0.015.png')
# ax_label(ax_a, 'a', x=0.9, color='white')
# show_png(ax_b, 'presentations/paper1/image_0.015_crop.png')
# ax_label(ax_b, 'b', x=0.9, color='white')
# show_png(ax_c, 'presentations/paper1/image_0.105.png')
# ax_label(ax_c, 'c', x=0.9, color='white')
# show_png(ax_d, 'presentations/paper1/image_0.105_crop.png')
# ax_label(ax_d, 'd', x=0.9, color='white')

# # fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
# # fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0, wspace=0)

# common.save_fig(fig, f'{path}/fig2.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig2.pdf', hide_metadata=True)


########################## fig 3 #########################
# fig = plt.figure(figsize=(w/2, 6))
# gs = matplotlib.gridspec.GridSpec(2, 1)#, height_ratios=(1, 1))
# gs0 = gs[0].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.4)
# gs1 = gs[1].subgridspec(1, 2, width_ratios=[1, 3], wspace=0.5)
# ax_a = fig.add_subplot(gs0[0])
# ax_b = fig.add_subplot(gs0[1])
# ax_c = fig.add_subplot(gs1[0])
# ax_d = fig.add_subplot(gs1[1])

# # gs = matplotlib.gridspec.GridSpec(2, 3, width_ratios=(0.3, 0.2, 0.5))
# # ax_a = fig.add_subplot(gs[0, :2])
# # ax_b = fig.add_subplot(gs[0, 2])
# # ax_c = fig.add_subplot(gs[1, 0])
# # ax_d = fig.add_subplot(gs[1, 1:])

# fig3_box_indices = [4, 11, 18, 25]
# box_counting.msd_combined.go(
#     ['eleanorlong001_no_overlap', 'eleanorlong001'],
#     ax_a,
#     box_size_indices=fig3_box_indices,
#     legend_num_boxes=True,
# )
# ax_a.set_ylim(0.4e-3, 0.3e3)
# ax_label(ax_a, 'a')

# D_of_L_mult_props = dict(
#     save_data=False,
#     title='',
#     show_nofit_cutoff=False,
#     box_size_indices=fig3_box_indices,
#     show_fits=False,
#     plot_C_N_squared=False,
#     rescale_x=None,
# )
# box_counting.D_of_L.go(
#     'eleanorlong001_no_overlap',
#     'nmsdfitinter',
#     ax=ax_b,
#     plot_color='tab:blue',
#     labels_on_plot=False, show_legend=False,
#     show_slope=True,
#     **D_of_L_mult_props,
# )
# box_counting.D_of_L.go(
#     'eleanorlong001',
#     'nmsdfitinter',
#     ax=ax_b,
#     plot_color='tab:orange',
#     labels_on_plot_font_color='black',
#     **D_of_L_mult_props,
# )
# ax_label(ax_b, 'b')

# # ax_label(ax_c, 'c', x=-0.3)
# ax_c.axis('off') 
# show_png2(fig, 'presentations/paper1/fig3_overlapped.png')

# box_counting.quantify_overlap_show.go(
#     'eleanorlong001',
#     ax_d,
# )
# ax_label(ax_d, 'd', x=0.3)

# common.save_fig(fig, f'{path}/fig3.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig3.pdf', hide_metadata=True)


# ########################## fig 4 #########################
# fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(w, 4))
# box_counting.D_of_L.go(
#     'eleanorlong010',
#     # 'var',
#     'nmsdfitinter',
#     ax=ax_a,
#     save_data=False,
#     title='',
#     show_nofit_cutoff=False,
#     max_num_boxes=6,
# )
# ax_a.set_xlim(0.5e-1, 1e6)
# # ax_a.yaxis.labelpad = -10
# # ax_a.xaxis.labelpad = -2
# ax_label(ax_a, 'a', x=0.95)

# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
# phi = 0.114
# D0_over_Sk = (1 + phi) / ( 1- phi)**3
# ax_b.hlines(D0_over_Sk, *ax_b.get_xlim(), color='gray', linestyle=':', label='$(1+\phi)/(1-\phi)^3$')
# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong001', 'eleanorlong010'],
#     ax=ax_b,
#     sources=['MSD_first', 'timescaleint_nmsdfitinter', 'timescaleint_nofit_cropped_nmsdfitinter'],
#     logarithmic_y=False,
#     legend_fontsize=7,
#     markers=[[None, 'o', 'x'], [None, 'o', 'x']]
# )
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_b.yaxis.set_major_locator(ticks_0p5)
# # ax_b.yaxis.labelpad = -12
# # ax_b.xaxis.labelpad = -2
# ax_label(ax_b, 'b')

# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong010', 'sim_nohydro_011_L320_longer_merged',
#     #  'brennan_hydro_010_L544'
#      ],
#     ax=ax_c,
#     sources=['MSD_first', 'D_of_L_theory', 'timescaleint_nmsdfitinter'],
#     logarithmic_y=False,
#     legend_fontsize=6,
#     discrete_colors=True,
#     markers=[[None, 'none', 'o'], [None, 'none', 'o'], [None, 'none', 'o']],
# )
# ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_c.yaxis.set_major_locator(ticks_0p5)
# ax_c.set_xlim(*DS_OVERLAPPED_XLIM)
# # ax_c.yaxis.labelpad = -12
# # ax_c.xaxis.labelpad = -2
# ax_label(ax_c, 'c')

# common.save_fig(fig, f'{path}/fig4.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig4.pdf', hide_metadata=True)


########################## fig 5 #########################
fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(w, 4))
box_counting.msd_single.go(
    'eleanorlong010',
    # show_timescaleint_replacement=True,
    # export_destination=f'{path}/fig6a.pdf',
    ax=ax_a,
    rescale_x=box_counting.msd_single.RESCALE_X_L2,
    rescale_y=box_counting.msd_single.RESCALE_Y_L2,
    legend_location='lower right',
    legend_fontsize=8,
    box_size_indices=[1, 6, 11, 16, 21, 25],
    show_nointer_theory_limits=True,
    # show_timescaleint_replacement=True,
)
# I tried linear y scale, it looks bad
ax_a.set_xlim(1e-2, 1e3)
ax_a.set_ylim(3e-3, 7e-2)

visualisation.Ds_overlapped_mult.go(
    ['eleanorlong001', 'eleanorlong010'],
    ax=ax_b,
    sources=['MSD_first', 'boxcounting_collective', 'D_of_L_theory'],
    logarithmic_y=False,
)
ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
ax_b.yaxis.set_major_locator(ticks_0p5)
ax_b.set_xlim(*DS_OVERLAPPED_XLIM)

visualisation.Ds_overlapped.go(
    'eleanorlong010',
    ax=ax_c,
    sources=['MSD_first', 'timescaleint_nmsdfitinter', 'boxcounting_collective', 'D_of_L_theory'],
    logarithmic_y=False,
)
ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
ax_c.yaxis.set_major_locator(ticks_0p5)
ax_c.set_xlim(*DS_OVERLAPPED_XLIM)
ax_label(ax_a, 'a')
ax_label(ax_b, 'b')
ax_label(ax_c, 'c')
common.save_fig(fig, f'{path}/fig5.png', hide_metadata=True)
common.save_fig(fig, f'{path}/fig5.pdf', hide_metadata=True)


########################## fig 6 #########################
# fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(w/2, 6))
# isf.show_Fs_overlayed.go(
#     'eleanorlong001',
#     ax_a,
#     target_ks = (0.2, 0.8, 2.4, 4.8),
#     SHOW_FIT=False
# )

# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong001'],
#     sources=['MSD_first', 'f_first_first', 'F_s_first'],
#     ax=ax_b,
#     logarithmic_y=False,
#     plot_against_k=True,
#     legend_fontsize=7,
# )
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_b.yaxis.set_major_locator(ticks_0p5)
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)

# # isf.show_Fs_overlayed.go(
# #     'eleanorlong001',
# #     ax_c,
# #     target_ks = [0.07],
# #     SHOW_FIT=False
# # )
# # ax_c.semilogy()
# # ax_c.set_xscale('linear')
# # ax_c.set_ylim(0.94, 1.01)
# # ax_c.set_xlim(-5, 150)
# # ax_c.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter()) # prevent scientific notation on axes
# # ax_c.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()) # prevent scientific notation on axes

# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong001_crop1.0', 'eleanorlong001_crop0.5', 'eleanorlong001_crop0.25', 'eleanorlong001_crop0.125',
#         #'eleanorlong001_crop0.0625'
#     ],
#     sources=['MSD_first', 'f_first_first'],
#     ax=ax_c,
#     plot_against_k=True,
#     legend_fontsize=7,
#     markers='o',
# )
# ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_c.yaxis.set_major_locator(ticks_0p5)
# ax_c.set_xlim(*DS_OVERLAPPED_XLIM)

# ax_label(ax_a, 'a', x=0.95)
# ax_label(ax_b, 'b', x=0.95)
# ax_label(ax_c, 'c', x=0.95)
# ax_label(ax_d, 'd', x=0.95)
# common.save_fig(fig, f'{path}/fig6.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig6.pdf', hide_metadata=True)


# ######################### fig 7 #########################
# fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(w/2, 3))
# DS_OVERLAPPED_YLIM_FKT = (0.8, 2.4)
# DS_OVERLAPPED_YLIM_FKT = DS_OVERLAPPED_YLIM

# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong001', 'eleanorlong010'],
#     ax=ax_a,
#     sources=['MSD_first', 'f_first_first'],
#     logarithmic_y=False,
#     legend_fontsize=7,
# )
# ax_a.set_ylim(*DS_OVERLAPPED_YLIM_FKT)
# ax_a.set_xlim(*DS_OVERLAPPED_XLIM)
# # ax_a.yaxis.labelpad = -5
# # ax_a.xaxis.labelpad = -2
# ax_label(ax_a, 'a')

# visualisation.Ds_overlapped_mult.go(
#     ['eleanorlong010', 'sim_nohydro_010_L640_div8', 'brennan_hydro_010_L544'],
#     ax=ax_b,
#     sources=['MSD_first', 'D0Sk_theory', 'f_first_first'],
#     logarithmic_y=False,
#     legend_fontsize=5,
#     theory_color='black',
#     discrete_colors=True,
#     file_labels=['experiment', 'sim. no hydro.', 'sim. with hydro.']
# )
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM_FKT)
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
# # ax_b.yaxis.labelpad = -5
# # ax_b.xaxis.labelpad = -2
# ax_label(ax_b, 'b')

# common.save_fig(fig, f'{path}/fig7.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig7.pdf', hide_metadata=True)


# ######################### fig 8 #########################
# fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(w/2, 3))

# visualisation.Ds_overlapped.go(
#     'eleanorlong010',
#     ax=ax_a,
#     sources=['MSD_first', 'f_first_first', 'f_long', 'timescaleint_nmsdfitinter'],
#     logarithmic_y=False,
#     legend_fontsize=6,
# )
# ax_a.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_a.yaxis.set_major_locator(ticks_0p5)
# ax_a.set_xlim(*DS_OVERLAPPED_XLIM)
# # ax_a.yaxis.labelpad = -5
# # ax_a.xaxis.labelpad = -2
# ax_label(ax_a, 'a')

# visualisation.Ds_overlapped.go(
#     'sim_nohydro_010_L640',
#     ax=ax_b,
#     sources=['MSD_first', 'f_first_first', 'f_short', 'f_long', 'timescaleint_nmsdfitinter'],
#     logarithmic_y=False,
#     legend_fontsize=6,
# )
# ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_b.yaxis.set_major_locator(ticks_0p5)
# ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
# # ax_b.yaxis.labelpad = -5
# # ax_b.xaxis.labelpad = -2
# ax_label(ax_b, 'b')

# common.save_fig(fig, f'{path}/fig8.png', hide_metadata=True)
# common.save_fig(fig, f'{path}/fig8.pdf', hide_metadata=True)