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
ONEROW_HEIGHT = 3.5
TWOROW_HEIGHT = 6.2

plt.rcParams.update({
    'axes.labelsize': 12,
    "legend.handlelength": 1.2, # space either side of coloured dot
    'legend.handletextpad': 0.4,
    'legend.handleheight': 0.7, # this is the default
    'legend.labelspacing': 0.4, # vertical spacing between labels
    'legend.borderpad': 0.4, # default
})

DS_OVERLAPPED_YLIM = (0.87, 2.5)
# DS_OVERLAPPED_YLIM_FKT = (0.8, 2.4)
DS_OVERLAPPED_YLIM_FKT = DS_OVERLAPPED_YLIM

DS_OVERLAPPED_XLIM = (0.05, 50)
# DS_OVERLAPPED_XLIM_K = (2*np.pi/50, 2*np.pi/0.05)
DS_OVERLAPPED_XLIM_K = (0.26, 45)
DS_OVERLAPPED_XLIM_K_FLIPPED = (2*np.pi/DS_OVERLAPPED_XLIM[1], 2*np.pi/DS_OVERLAPPED_XLIM[0])

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

LINESTYLE_FKT = (0, (1, 0.8))
LINESTYLE_COUNTING = 'dashed'

SIM_EXP_TITLE_FONTSIZE = 15

PHI_FILE_LABELS = ['$\phi=0.02$', '$\phi=0.11$']

if __name__ == '__main__':
    pass
    
    ######################### fig nmsd_fit #########################
    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(WIDTH/2, TWOROW_HEIGHT),
                                                                    # constrained_layout=True,
                                                    #   gridspec_kw={'wspace': 0.4, 'hspace': 0.5},
                                                    sharey='row'
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
        box_size_indices=[7, 11, 15, 19],
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
    
    mult_exponent = -3
    ax_a.text(-0.24, 1.06, rf'$\times 10^{{{mult_exponent}}}$', transform=ax_a.transAxes)
    ax_a.yaxis.labelpad = -5

    wanted_ticks = [0.6e-2, 1e-2, 4e-2]
    def tick_formatter(x, pos):
        assert False
        print('doin it')
        if x in wanted_ticks:
            mult = 10**mult_exponent
            return f'{x*mult}'
        else:
            return ''

    box_counting.msd_single.go(
        'eleanorlong010',
        ax=ax_b,
        legend_location='lower right',
        disable_ylabel=True,
        **msd_zoom_props,
    )
    
    
    ax_a.yaxis.set_major_formatter(tick_formatter)
    ax_a.yaxis.set_minor_formatter(tick_formatter)
    
    # ax_b.get_legend().remove()
    # ax_b.set_xlim(*MSD_ZOOM_XLIM)
    # ax_b.set_ylim(*MSD_ZOOM_YLIM)
    # # ax_b.set_yscale('linear') # I tried linear y scale, it looks bad
    # # do_ticks(ax_b)

    # Ds_mult_props = dict(
    #     logarithmic_y=False,
    #     legend_fontsize=LEGEND_FONTSIZE,
    #     file_labels=PHI_FILE_LABELS,
    #     source_labels=['Eq. ? fit', 'theory'],
    #     colors=[[COLOR_PHI002, COLOR_PHI002_THEORY], [COLOR_PHI011, COLOR_PHI011_THEORY]],
    #     sources=[
    #         'boxcounting_collective_var',
    #         # 'boxcounting_collective_nmsdfitinter',
    #         'D_of_L_theory'
    #     ],
    #     linestyles=[['none', LINESTYLE_COUNTING]]*2,
    # )

    # visualisation.Ds_overlapped_mult.go(
    #     ['sim_nohydro_002_L320_longer_merged', 'sim_nohydro_011_L320_longer_merged'],
    #     ax=ax_c,
    #     **Ds_mult_props
    # )
    # ax_c.set_ylim(*DS_OVERLAPPED_YLIM)
    # ax_c.yaxis.set_major_locator(ticks_0p5)
    # ax_c.set_xlim(*DS_OVERLAPPED_XLIM)

    # visualisation.Ds_overlapped_mult.go(
    #     ['eleanorlong001', 'eleanorlong010'],
    #     ax=ax_d,
    #     disable_ylabel=True,
    #     **Ds_mult_props
    # )
    # ax_d.set_ylim(*DS_OVERLAPPED_YLIM)
    # ax_d.yaxis.set_major_locator(ticks_0p5)
    # ax_d.set_xlim(*DS_OVERLAPPED_XLIM)

    # ax_label(ax_b, 'b')
    # ax_label(ax_c, 'c', x=0.95)
    # ax_label(ax_d, 'd', x=0.95)
    common.save_fig(fig, f'{path}/fig_nmsd_fit.png', hide_metadata=True)
    common.save_fig(fig, f'{path}/fig_nmsd_fit.pdf', hide_metadata=True)

