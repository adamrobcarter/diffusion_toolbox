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
import isf.show_S_of_k

import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec
import matplotlib.colors
import numpy as np
import cmocean

"""
Figures should be prepared with the PDF layout in mind.
Individual figures should not be longer than one page and with a width that corresponds to
1 column (85 mm) or 2 columns (180 mm).

All images must have a resolution of 300 dpi at final size.
"""

path = 'presentations/paper1'
WIDTH = 12
ONEROW_HEIGHT = 3.1
TWOROW_HEIGHT = 6.2

plt.rcParams.update({
    'axes.labelsize': 12,
    "legend.handlelength": 1.2, # space either side of coloured dot
    'legend.handletextpad': 0.4,
    'legend.handleheight': 0.7, # this is the default
    'legend.labelspacing': 0.4, # vertical spacing between labels
    'legend.borderpad': 0.4, # default
    'legend.framealpha': 0.96,
})

DS_OVERLAPPED_YLIM = (0.87, 2.5)
# DS_OVERLAPPED_YLIM_FKT = (0.8, 2.4)
DS_OVERLAPPED_YLIM_FKT = DS_OVERLAPPED_YLIM
DS_OVERLAPPED_YLIM_FKT_SMALL = (0.9, 2)

DS_OVERLAPPED_XLIM = (0.05, 50)
# DS_OVERLAPPED_XLIM_K = (2*np.pi/50, 2*np.pi/0.05)
DS_OVERLAPPED_XLIM_K = (0.26, 45)
DS_OVERLAPPED_XLIM_L_FOR_FKT = (0.13, 50)

ticks_0p5 = plt.MultipleLocator(base=0.5)

LEGEND_FONTSIZE=8

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

def show_png(ax, file, border=None):
    image = plt.imread(file)
    ax.imshow(image, interpolation='none')
    if border:
        ax.axes.xaxis.set_visible(False) # you can use these two to hide the axes
        ax.axes.yaxis.set_visible(False) # but keep a black border around the plot
        for spine in ax.spines.values():
            spine.set_edgecolor(border)
            spine.set_linewidth(3)
    else:
        ax.axis('off') 

def show_png_offaxis(fig, file, where):
    image = plt.imread(file)
    print(image.shape)
    # fig.figimage(image, xo=50, yo=50,  origin='lower')

    ax = fig.add_axes(where)#, frameon=True)
    ax.imshow(image, interpolation='none')
    # ax.imshow(image, interpolation='none')
    # ax.axis('off') 
    ax.axes.xaxis.set_visible(False) # you can use these two to hide the axes
    ax.axes.yaxis.set_visible(False) # but keep a black border around the plot

    return ax

def set_title_sim(ax):
    ax.set_title('simulation', fontsize=SIM_EXP_TITLE_FONTSIZE, pad=10)

def set_title_exp(ax):
    ax.set_title('experiment', fontsize=SIM_EXP_TITLE_FONTSIZE, pad=10)

def truncate_colormap(cmap, minval=0.0, maxval=1.0): # https://stackoverflow.com/a/18926541/1103752
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, 200)))
    return new_cmap

PHI002_CMAP_POINT = 0.3
COLOR_PHI002 = cmocean.cm.ice(PHI002_CMAP_POINT)
COLOR_PHI011 = 'tab:pink'
COLOR_PHI002_THEORY = '#202346'
COLOR_PHI011_THEORY = '#551140' 
# COLOR_NONOVERLAPPED = 'tab:green'
COLOR_NONOVERLAPPED = 'limegreen'
COLORMAP_FKT_CROPS = cmocean.tools.crop_by_percent(cmocean.tools.crop_by_percent(cmocean.cm.ice, PHI002_CMAP_POINT*100, which='min', N=None), 25, which='max', N=None)
# COLORMAP_BOXES = cmocean.cm.solar
# COLORMAP_BOXES = matplotlib.cm.inferno
# COLORMAP_BOXES = truncate_colormap(matplotlib.cm.afmhot, 0.2, 0.77)
COLORMAP_BOXES = common.colormap
# COLORMAP_KS = lambda a : common.colormap(0.75-a) # 0.75 is a hack cause there are 4 plots on the graph (n-1)/n
COLORMAP_KS = common.colormap

MARKER_FKT = 'd'
MARKER_COUNTING = 'o'
MARKER_PHENOMFIT = 's'

LINESTYLE_FKT = (0, (1, 0.8))
LINESTYLE_COUNTING = 'dashed'

SIM_EXP_TITLE_FONTSIZE = 15

PHI_FILE_LABELS = ['$\phi=0.02$', '$\phi=0.11$']

GRIDSPEC_WSPACE = 0.07
GRIDSPEC_HSPACE = 0.27

if __name__ == '__main__':
    pass
    ########################## fig 1 #########################
    fig = plt.figure(figsize=(WIDTH/2, 5))
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
    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(WIDTH/2, 5.15),
                                                     constrained_layout=True,
                                                     gridspec_kw={'wspace': 0.09, 'hspace': 0})

    show_png(ax_a, 'presentations/paper1/image_0.015.png', border=COLOR_PHI002)
    ax_label(ax_a, 'a', x=0.9, color='white')
    show_png(ax_b, 'presentations/paper1/image_0.015_crop.png', border=COLOR_PHI002)
    ax_label(ax_b, 'b', x=0.9, color='white')
    show_png(ax_c, 'presentations/paper1/image_0.105.png', border=COLOR_PHI011)
    ax_label(ax_c, 'c', x=0.9, color='white')
    show_png(ax_d, 'presentations/paper1/image_0.105_crop.png', border=COLOR_PHI011)
    ax_label(ax_d, 'd', x=0.9, color='white')

    # fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    # fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0, wspace=0)

    common.save_fig(fig, f'{path}/fig2.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig2.pdf', hide_metadata=True)


    ######################### fig 3 #########################
    fig = plt.figure(figsize=(WIDTH/2, 6))
    gs = matplotlib.gridspec.GridSpec(2, 1)#, height_ratios=(1, 1))
    gs0 = gs[0].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.4)
    gs1 = gs[1].subgridspec(1, 2, width_ratios=[1, 3], wspace=0.5)
    ax_a = fig.add_subplot(gs0[0])
    ax_b = fig.add_subplot(gs0[1])
    ax_c = fig.add_subplot(gs1[0])
    ax_d = fig.add_subplot(gs1[1])

    # gs = matplotlib.gridspec.GridSpec(2, 3, width_ratios=(0.3, 0.2, 0.5))
    # ax_a = fig.add_subplot(gs[0, :2])
    # ax_b = fig.add_subplot(gs[0, 2])
    # ax_c = fig.add_subplot(gs[1, 0])
    # ax_d = fig.add_subplot(gs[1, 1:])

    fig3_box_indices = [4, 11, 18, 25]
    box_counting.msd_combined.go(
        ['eleanorlong001_no_overlap', 'eleanorlong001'],
        ax_a,
        box_size_indices=fig3_box_indices,
        legend_num_boxes=True,
        colors=[COLOR_NONOVERLAPPED, COLOR_PHI002]
    )
    ax_a.set_ylim(0.4e-3, 0.3e3)
    ax_label(ax_a, 'a')
    ax_a.yaxis.labelpad = -1
    # ax_a.legend(loc='upper center', fontsize=LEGEND_FONTSIZE)

    D_of_L_mult_props = dict(
        save_data=False,
        title='',
        show_nofit_cutoff=False,
        box_size_indices=fig3_box_indices,
        show_short_fits=False,
        show_long_fits=False,
        plot_C_N_squared=False,
        rescale_x=None,
    )
    box_counting.D_of_L.go(
        'eleanorlong001_no_overlap',
        'var',# 'nmsdfitinter',
        ax=ax_b,
        plot_color=COLOR_NONOVERLAPPED,
        labels_on_plot=False, show_legend=False,
        show_slope=True,
        **D_of_L_mult_props,
    )
    box_counting.D_of_L.go(
        'eleanorlong001',
        'var',# 'nmsdfitinter',
        ax=ax_b,
        plot_color=COLOR_PHI002,
        labels_on_plot_font_color='black',
        **D_of_L_mult_props,
    )
    ax_label(ax_b, 'b')
    ax_b.yaxis.labelpad = -1

    # ax_label(ax_c, 'c', x=-0.3)
    ax_c.axis('off')
    y_offset = 0.05
    x_offset = 0.02
    width  = 0.33
    height = 0.22
    ax = show_png_offaxis(fig, 'presentations/paper1/fig3_nonOver.png', [x_offset, y_offset+height, width, height])
    _  = show_png_offaxis(fig, 'presentations/paper1/fig3_Over.png',    [x_offset, y_offset,        width, height])
    
    ax_label(ax, 'c')

    box_counting.quantify_overlap_show.go(
        'eleanorlong001',
        ax_d,
    )
    ax_label(ax_d, 'd', x=0.3)
    ax_d.legend(fontsize=LEGEND_FONTSIZE)

    common.save_fig(fig, f'{path}/fig3.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig3.pdf', hide_metadata=True)


    ########################## fig tsi #########################
    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2,
                                                     figsize=(WIDTH/2, TWOROW_HEIGHT),
                                                     sharey='row',
                                                     gridspec_kw={'wspace': GRIDSPEC_WSPACE, 'hspace': GRIDSPEC_HSPACE},
                                                     )
    set_title_sim(ax_a)
    set_title_exp(ax_b)

    D_of_L_props = dict(
        save_data=False,
        title='',
        show_nofit_cutoff=False,
        max_num_boxes=6,
        colormap=COLORMAP_BOXES,
        show_short_fits=True,
        show_long_fits=False,
    )

    box_counting.D_of_L.go(
        'sim_nohydro_011_L640_longer_merged',
        'var', #'nmsdfitinter',
        ax=ax_a,
        **D_of_L_props,
    )
    ax_a.set_xlim(0.5e-1, 1e6)
    ax_a.yaxis.labelpad = -1
    ax_a.xaxis.labelpad = -1
    ax_label(ax_a, 'a', x=0.95)

    box_counting.D_of_L.go(
        'eleanorlong010',
        'var', # 'nmsdfitinter',
        ax=ax_b,
        disable_ylabel=True,
        **D_of_L_props,
    )
    ax_b.set_xlim(0.5e-1, 1e6)
    ax_b.xaxis.labelpad = -1
    ax_label(ax_b, 'b', x=0.95)

    # ax_c.set_xlim(*DS_OVERLAPPED_XLIM)
    phi = 0.114**3
    D0_over_Sk = (1 + phi) / ( 1- phi)
    # ax_c.hlines(D0_over_Sk, *ax_c.get_xlim(), color='gray', linestyle=':', label='$(1+\phi)/(1-\phi)^3$')

    Ds_ov_mult_props = dict(
        logarithmic_y=False,
        legend_fontsize=LEGEND_FONTSIZE-1,
        file_labels=PHI_FILE_LABELS,
        source_labels=['theory', 'tsi no fit'],
        colors=[[COLOR_PHI002_THEORY, COLOR_PHI002], [COLOR_PHI011_THEORY, COLOR_PHI011]],
        markers=[['none', MARKER_COUNTING]] * 2,
        linestyles=[[LINESTYLE_COUNTING, 'none']]*2,
    )

    visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L640_longer_merged', 'sim_nohydro_011_L640_longer_mergedD'],
        ax=ax_c,
        # sources=['D_of_L_theory', 'timescaleint_nmsdfitinter', 'timescaleint_nofit_cropped_nmsdfitinter'],
        sources=['D_of_L_theory', 'timescaleint_nofit_cropped_var'],
        **Ds_ov_mult_props
        # markers=[[None, MARKER_COUNTING, 'x'], [None, MARKER_COUNTING, 'x']]
    )
    ax_c.set_xlim(*DS_OVERLAPPED_XLIM)
    ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_c.yaxis.set_major_locator(ticks_0p5)
    # ax_b.yaxis.labelpad = -12
    ax_c.xaxis.labelpad = -1
    ax_label(ax_c, 'c')#, x=0.95)
    ax_c.get_legend().remove()

    visualisation.Ds_overlapped_mult.go(
        ['eleanorlong001', 'eleanorlong010',
        #  'brennan_hydro_010_L544'
         ],
        ax=ax_d,
        # sources=['D_of_L_theory', 'timescaleint_nmsdfitinter', 'timescaleint_nofit_cropped_nmsdfitinter'],
        sources=['D_of_L_theory', 'timescaleint_nofit_cropped_var'],
        **Ds_ov_mult_props,
        disable_ylabel=True,
        # markers=[[None, 'none', 'o'], [None, 'none', 'o'], [None, 'none', 'o']],
    )
    ax_d.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_d.yaxis.set_major_locator(ticks_0p5)
    ax_d.set_xlim(*DS_OVERLAPPED_XLIM)
    ax_d.xaxis.labelpad = -1
    ax_label(ax_d, 'd', x=0.95)
    ax_d.get_legend().set_bbox_to_anchor((-0.13, 0, 1, 1))

    common.save_fig(fig, f'{path}/fig_tsi.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig_tsi.pdf', hide_metadata=True)


    ######################### fig nmsd_fit #########################
    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(WIDTH/2, TWOROW_HEIGHT),
                                                                    # constrained_layout=True,
                                                    #   gridspec_kw={'wspace': 0.4, 'hspace': 0.5},
                                                    sharey='row',
                                                    gridspec_kw={'wspace': GRIDSPEC_WSPACE, 'hspace': GRIDSPEC_HSPACE},
                                                    )
    
    set_title_sim(ax_a)
    set_title_exp(ax_b)
    MSD_ZOOM_XLIM = (10e-2, 6e2)
    MSD_ZOOM_YLIM = (5e-3, 4e-2)
    msd_zoom_props = dict(
        rescale_x=box_counting.msd_single.RESCALE_X_L2,
        rescale_y=box_counting.msd_single.RESCALE_Y_L2,
        legend_fontsize=LEGEND_FONTSIZE,
        # box_size_indices=[1, 6, 11, 16, 21, 25],
        # box_size_indices=[3, 7, 11, 15, 19, 23],
        box_size_indices=[3, 7, 11, 15, 19, 23],
        show_nointer_theory_limits=True,
        nointer_theory_limit_labels=['$D_\mathrm{self}$ limit', '$D_\mathrm{coll}$ limit']
    )

    box_counting.msd_single.go(
        'sim_nohydro_011_L320_longer_merged',
        ax=ax_a,
        legend_location='lower right',
        **msd_zoom_props
    )
    ax_a.set_xlim(*MSD_ZOOM_XLIM)
    ax_a.set_ylim(*MSD_ZOOM_YLIM)
    ax_label(ax_a, 'a')
    # ax_a.set_yscale('linear') # I tried linear y scale, it looks bad
    ax_a.yaxis.labelpad = -1
    ax_a.xaxis.labelpad = -1
    # ax_a.get_legend().set_bbox_to_anchor((0, 0, 1.13, 1))
    ax_a.get_legend().remove()
    
    def do_ticks(ax):
        mult_exponent = -2
        ax.text(-0.23, 1.06, rf'$\times 10^{{{mult_exponent}}}$', transform=ax.transAxes)
        # ax.yaxis.labelpad = -5

        wanted_ticks = [0.5e-2, 1e-2, 2e-2, 4e-2]
        def tick_formatter(x, pos):
            if x in wanted_ticks:
                mult = 10**mult_exponent
                return f'{x/mult:.1f}'
            else:
                return ''
        
        ax.yaxis.set_major_formatter(tick_formatter)
        ax.yaxis.set_minor_formatter(tick_formatter)

    box_counting.msd_single.go(
        'eleanorlong010',
        ax=ax_b,
        legend_location='lower right',
        disable_ylabel=True,
        show_second_legend=False,
        **msd_zoom_props,
    )
    # ax_b.get_legend().remove()
    ax_b.set_xlim(*MSD_ZOOM_XLIM)
    ax_b.set_ylim(*MSD_ZOOM_YLIM)
    # ax_b.set_yscale('linear') # I tried linear y scale, it looks bad
    # do_ticks(ax_b)
    ax_b.xaxis.labelpad = -1
    ax_label(ax_b, 'b')
    ax_b.get_legend().set_bbox_to_anchor((-0.95, 0, 1, 1))

    Ds_mult_props = dict(
        logarithmic_y=False,
        legend_fontsize=LEGEND_FONTSIZE,
        file_labels=PHI_FILE_LABELS,
        source_labels=['Eq. 6 fit', 'theory'],
        colors=[[COLOR_PHI002, COLOR_PHI002_THEORY], [COLOR_PHI011, COLOR_PHI011_THEORY]],
        sources=[
            'boxcounting_collective_var',
            # 'boxcounting_collective_nmsdfitinter',
            'D_of_L_theory'
        ],
        linestyles=[['none', LINESTYLE_COUNTING]]*2,
        markers = MARKER_PHENOMFIT
    )

    
    # matplotlib is being weird and this has to be after we do ax_b
    # probably something to do with sharey
    do_ticks(ax_a)

    D0_over_Sk0 = (1 + 0.114) / ( 1 - 0.114)**3
    D0_over_Sk0_2 = (1 + 0.016) / ( 1 - 0.016)**3
    hlines_kwargs = dict(
        linestyle='dotted',
        label='$D_\mathrm{coll}$',
        color=COLOR_PHI011_THEORY,
        linewidth=1,
    )
    ax_c.hlines(D0_over_Sk0, *DS_OVERLAPPED_XLIM, **hlines_kwargs)
    visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L320', 'sim_nohydro_011_L320'],
        ax=ax_c,
        **Ds_mult_props
    )
    ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_c.yaxis.set_major_locator(ticks_0p5)
    ax_c.set_xlim(*DS_OVERLAPPED_XLIM)
    ax_c.yaxis.labelpad = -1
    ax_c.xaxis.labelpad = -1
    ax_c.get_legend().remove()
    ax_label(ax_c, 'c')
    

    ax_d.hlines(D0_over_Sk0, *DS_OVERLAPPED_XLIM, **hlines_kwargs)
    visualisation.Ds_overlapped_mult.go(
        ['eleanorlong001', 'eleanorlong010'],
        ax=ax_d,
        disable_ylabel=True,
        **Ds_mult_props
    )
    ax_d.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_d.yaxis.set_major_locator(ticks_0p5)
    ax_d.set_xlim(*DS_OVERLAPPED_XLIM)
    ax_d.xaxis.labelpad = -1
    ax_d.get_legend().set_bbox_to_anchor((-0.12, 0, 1, 1))
    ax_label(ax_d, 'd', x=0.95)
    
    common.save_fig(fig, f'{path}/fig_nmsd_fit.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig_nmsd_fit.pdf', hide_metadata=True)


    ######################### fig 6 (fkt) #########################
    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2,
                                                     figsize=(WIDTH/2, 6),
                                                     gridspec_kw={'wspace': 0.28, 'hspace': 0.27},
                                                     )
    isf.show_Fs_overlayed.go(
        'eleanorlong001',
        ax_a,
        target_ks = (0.2, 0.8, 2.4, 4.8),
        SHOW_FIT=False,
        colormap=COLORMAP_KS
    )
    ax_a.yaxis.labelpad = -1
    ax_a.xaxis.labelpad = -1

    visualisation.Ds_overlapped_mult.go(
        ['eleanorlong001'],
        sources=['f_first_first', 'F_s_first'],
        ax=ax_b,
        logarithmic_y=False,
        plot_against_k=True,
        legend_fontsize=LEGEND_FONTSIZE,
        source_labels=['$f(k, t)$', '$F_s(k, t)$'],
        file_labels=[' '],
        colors=[[COLOR_PHI002, COLOR_PHI002]],
        markers=[[MARKER_FKT, 'x']]
    )
    ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_b.yaxis.set_major_locator(ticks_0p5)
    ax_b.set_xlim(*DS_OVERLAPPED_XLIM_K)
    ax_b.yaxis.labelpad = -1
    ax_b.xaxis.labelpad = -1

    # isf.show_Fs_overlayed.go(
    #     'eleanorlong001',
    #     ax_c,
    #     target_ks = [0.07],
    #     SHOW_FIT=False
    # )
    # ax_c.semilogy()
    # ax_c.set_xscale('linear')
    # ax_c.set_ylim(0.94, 1.01)
    # ax_c.set_xlim(-5, 150)
    # ax_c.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter()) # prevent scientific notation on axes
    # ax_c.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()) # prevent scientific notation on axes
    L = np.array([361.6])
    get_L = lambda c : rf'$L_x={L[0]*c:.0f}\mathrm{{\mu m}}$'
    visualisation.Ds_overlapped_mult.go(
        ['eleanorlong001_crop1.0', 'eleanorlong001_crop0.5', 'eleanorlong001_crop0.25', 'eleanorlong001_crop0.125',
            #'eleanorlong001_crop0.0625'
        ],
        sources=['f_first_first'],
        ax=ax_c,
        plot_against_k=True,
        legend_fontsize=LEGEND_FONTSIZE,
        markers='d',
        file_labels=[get_L(1), get_L(0.5), get_L(0.25), get_L(0.125)],
        source_labels=['', ''],
        colors=[[COLORMAP_FKT_CROPS(0.0)], [COLORMAP_FKT_CROPS(0.33)], [COLORMAP_FKT_CROPS(0.66)], [COLORMAP_FKT_CROPS(1.0)]]
    )
    ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_c.yaxis.set_major_locator(ticks_0p5)
    ax_c.set_xlim(*DS_OVERLAPPED_XLIM_K)
    ax_c.xaxis.labelpad = -5
    ax_c.yaxis.labelpad = -1

    visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L320', 'sim_nohydro_002_L640_crop320'],
        sources=[
            'f_first_first'
        ],
        ax=ax_d,
        plot_against_k=True,
        markers=MARKER_FKT,
        file_labels=['periodic boundary', 'non periodic boundary'],
        source_labels=[''],
        colors=[['olivedrab'], ['darkgreen']]
    )
    ax_d.text(0.95, 0.85, 'simulation', ha='right', fontsize=15, transform=ax_d.transAxes)
    ax_d.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_d.yaxis.set_major_locator(ticks_0p5)
    ax_d.set_xlim(*DS_OVERLAPPED_XLIM_K)
    ax_d.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1, 0.8)) # puts upper right at specified position
    ax_d.xaxis.labelpad = -5
    ax_d.yaxis.labelpad = -1

    ax_label(ax_a, 'a', x=0.95)
    ax_label(ax_b, 'b')#, x=0.95)
    ax_label(ax_c, 'c')#, x=0.95)
    ax_label(ax_d, 'd')#, x=0.95)
    common.save_fig(fig, f'{path}/fig6.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig6.pdf', hide_metadata=True)


    ######################## fig 7 #########################
    fig, (ax_a, ax_b) = plt.subplots(1, 2,
                                     figsize=(WIDTH/2, ONEROW_HEIGHT),
                                     sharey='row',
                                     gridspec_kw={'wspace': GRIDSPEC_WSPACE},
                                     )
    set_title_sim(ax_a)
    set_title_exp(ax_b)

    fig7_props = dict(
        logarithmic_y=False,
        legend_fontsize=LEGEND_FONTSIZE,
        colors=[[COLOR_PHI002, COLOR_PHI002_THEORY], [COLOR_PHI011, COLOR_PHI011_THEORY]],
        file_labels=PHI_FILE_LABELS,
        sources=['f_first_first', 'D0Sk_theory'],
        markers=[[MARKER_FKT, 'none']]*2,
        linestyles=[['none', LINESTYLE_FKT]]*2,
        source_labels=['$f(k, t)$', 'theory']
    )

    visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_002_L320_longer_mergedD', 'sim_nohydro_011_L320_longer_mergedD'],#, 'brennan_hydro_011_L320'],
        #                                          ^^^^ atm brennnan does not need merging cause dt16 is no longer than dt0.5
        ax=ax_a,
        **fig7_props
    )
    ax_a.set_ylim(*DS_OVERLAPPED_YLIM_FKT_SMALL)
    ax_a.set_xlim(*DS_OVERLAPPED_XLIM_L_FOR_FKT)
    # ax_a.yaxis.labelpad = -5
    ax_a.xaxis.labelpad = -5
    ax_label(ax_a, 'a')#, x=0.85)
    ax_a.get_legend().remove()

    visualisation.Ds_overlapped_mult.go(
        ['eleanorlong001', 'eleanorlong010'],
        ax=ax_b,
        disable_ylabel=True,
        **fig7_props
    )
    ax_b.set_ylim(*DS_OVERLAPPED_YLIM_FKT_SMALL)
    ax_b.set_xlim(*DS_OVERLAPPED_XLIM_L_FOR_FKT)
    # ax_b.yaxis.labelpad = -5
    ax_b.xaxis.labelpad = -5
    ax_label(ax_b, 'b', x=0.95)
    ax_b.get_legend().set_bbox_to_anchor((-0.13, 0, 1, 1))

    common.save_fig(fig, f'{path}/fig7.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig7.pdf', hide_metadata=True)


    ######################### fig 8 #########################
    fig, (ax_a, ax_b) = plt.subplots(1, 2,
                                     figsize=(WIDTH/2, ONEROW_HEIGHT),
                                     sharey='row',
                                     gridspec_kw={'wspace': GRIDSPEC_WSPACE},
                                     )
    set_title_sim(ax_a)
    set_title_exp(ax_b)

    fig8_props = dict(
        logarithmic_y=False,
        legend_fontsize=LEGEND_FONTSIZE,
        file_labels=[''],
        # sources=      ['timescaleint_nmsdfitinter', 'D_of_L_theory',      'f_first_first', 'D0Sk_theory'],
        sources=      ['timescaleint_nofit_cropped_var',          'D_of_L_theory',      'f_first_first', 'D0Sk_theory'],
        colors=       [[COLOR_PHI011,               COLOR_PHI011_THEORY,  COLOR_PHI011,    COLOR_PHI011_THEORY]],
        markers=      [[MARKER_COUNTING,            'none',               MARKER_FKT,      'none', ]],
        linestyles=   [['none',                     LINESTYLE_COUNTING,   'none',          LINESTYLE_FKT]],
        source_labels=['Countoscope',               'Countoscope theory', '$f(k, t)$',     '$f(k, t)$ theory']
    )

    visualisation.Ds_overlapped_mult.go(
        ['sim_nohydro_011_L320_longer_mergedD'],
        ax=ax_a,
        **fig8_props
    )
    ax_a.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_a.yaxis.set_major_locator(ticks_0p5)
    ax_a.set_xlim(*DS_OVERLAPPED_XLIM)
    ax_a.yaxis.labelpad = -1
    ax_a.xaxis.labelpad = -1
    ax_label(ax_a, 'a')#, x=0.95)
    ax_a.get_legend().remove()

    visualisation.Ds_overlapped_mult.go(
        ['eleanorlong010'],
        ax=ax_b,
        disable_ylabel=True,
        **fig8_props
    )
    ax_b.set_ylim(*DS_OVERLAPPED_YLIM)
    ax_b.yaxis.set_major_locator(ticks_0p5)
    ax_b.set_xlim(*DS_OVERLAPPED_XLIM)
    ax_b.xaxis.labelpad = -1
    ax_label(ax_b, 'b', x=0.95)
    ax_b.get_legend().set_bbox_to_anchor((-0.4, 0, 1, 1))

    common.save_fig(fig, f'{path}/fig8.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig8.pdf', hide_metadata=True)

    ############### S(k) appendix ###############
    fig, ax = plt.subplots(1, 1,
                                     figsize=(WIDTH/3, ONEROW_HEIGHT),
                                     )

    isf.show_S_of_k.go(
        'eleanorlong010',
        ax
    )

    common.save_fig(fig, f'{path}/fig_Sk_appendix.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig_Sk_appendix.pdf', hide_metadata=True)