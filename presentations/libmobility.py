import common
import visualisation.Ds_overlapped_mult
import isf.show_S_of_k
import matplotlib.pyplot as plt
import numpy as np
import isf.show_Fs_overlayed
import matplotlib.image
import matplotlib.offsetbox
import matplotlib.cm
import matplotlib.colors

figs = 'presentations/libmobility'

WALL_LONG_COLOR = 'royalblue'
WALL_SHORT_COLOR = 'cornflowerblue'
WALL_LONG_COLOR = matplotlib.cm.Blues(0.8)
WALL_SHORT_COLOR = matplotlib.cm.Blues(0.64)
OPEN_LONG_COLOR = 'indianred'
OPEN_SHORT_COLOR = 'lightcoral'

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern"
})

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

from presentations.presentations import show_png_on_axis


def show_svg_on_axis(ax, file, coords, size=1.0, color='black', data_coords=True):
    """
    Display an SVG file on a Matplotlib axis as a vector graphic (not rasterized).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    file : str
        Path to the .svg file.
    coords : tuple
        (x, y) coordinates where the SVG will be placed.
    size : float
        Scaling factor for the SVG.
    color : str
        Border color of the annotation box.
    data_coords : bool
        If True, coords are in data space; otherwise in axes fraction.
    """
    # Load the SVG as a vector object
    svg_box = matplotlib.offsetbox.SVG(file, zoom=size)

    ab = matplotlib.offsetbox.AnnotationBbox(
        svg_box,
        coords,
        xycoords="data" if data_coords else "axes fraction",
        frameon=True
    )

    # Customize the frame (optional)
    frame = ab.patch
    frame.set_edgecolor(color)

    ax.add_artist(ab)

fig, axs = plt.subplot_mosaic(
    """
    ac
    bc
    """,
    gridspec_kw=dict(width_ratios=(1, 2.2)),
    figsize=(6, 3.8)
)

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
        dict(
            file='ld_hydro_nbody_0.114_L2560_t450_1',
            source='D0Sk_theory',
            label='theory, no hydrodynamics',
            color='black'
        ),
    ),
    ax = axs['c'],
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
    allow_rescale_x=True,
)
# ax.set_ylim(0.5, 4)
axs['c'].set_ylim(0.8, 9)
axs['c'].set_xlim(6.3e-2, 3e1)
axs['c'].semilogy()
common.add_exponential_index_indicator(axs['c'], -1, (4.5e-1, 3.6e0), 'k', x_limits=(1e-2, 6e-1), )
# axs['c'].grid()
ax_label(axs['c'], 'c', x=0.95, y=0.98)
show_png_on_axis(axs['c'], 'presentations/libmobility/side_view_2.png',         (0.85, 0.5 ), 0.04, color=WALL_LONG_COLOR, data_coords=False)
show_png_on_axis(axs['c'], 'presentations/libmobility/side_view_potential.png', (0.85, 0.67), 0.04, color=OPEN_LONG_COLOR, data_coords=False)
# show_svg_on_axis(axs['c'], 'presentations/libmobility/side_view-one_wall.svg',         (0.85, 0.5 ), 0.04, color=WALL_LONG_COLOR, data_coords=False)
# show_svg_on_axis(axs['c'], 'presentations/libmobility/side_view-potential.png', (0.85, 0.67), 0.04, color=OPEN_LONG_COLOR, data_coords=False)
axs['c'].yaxis.labelpad = -8
axs['c'].xaxis.labelpad = -2


isf.show_S_of_k.go(
    'ld_hydro_nbody_0.114_L2560',
    axs['a'],
    source='F',
    show_realspace_axis=False,
    data_color=WALL_LONG_COLOR,
)
axs['a'].set_ylim(0.55, 1.15)
ax_label(axs['a'], 'a')
axs['a'].xaxis.labelpad = -2

isf.show_Fs_overlayed.go(
    'ld_hydro_nbody_0.114_L2560_t1h_1_fulltime',
    file_i=0,
    ax=axs['b'],
    target_ks=np.logspace(np.log10(1), np.log10(0.08), num=3),
    colormap=matplotlib.colors.LinearSegmentedColormap.from_list(
        "trunc_viridis", matplotlib.cm.Blues(np.linspace(0.3, 1, 256))
    ),
    show_fit=True,
)
ax_label(axs['b'], 'b', x=0.95)
axs['b'].xaxis.labelpad = -5


common.save_fig(fig, f'{figs}/collective_diffusion.pdf', hide_metadata=True)
common.save_fig(fig, f'{figs}/collective_diffusion.png', hide_metadata=True, dpi=300)

#### big dataset
fig, ax = plt.subplots()

msd_file = 'ld_hydro_dpstokes_0.114_L1280_short_unwrap_2d'
visualisation.Ds_overlapped_mult.go(
    (
        # short one
        dict(
            file='ld_hydro_dpstokes_0.114_L1280_short',
            source='f_t1',
            plot_index=np.index_exp[23:],
            color=WALL_SHORT_COLOR,
            msd_file=msd_file,
            label='single wall, $t=1\mathrm{{s}}$'
        ),
        # long one
        dict(
            file='ld_hydro_dpstokes_0.114_L1280',
            source='f_t1024',
            msd_file=msd_file,
            color=WALL_LONG_COLOR,
            label='single wall, $t=1024\mathrm{{s}}$'
        ),
        # dict(
        #     file='ld_dpstokes_0.114_nowalls_L1280_t8h_32s',
        #     source='f_t64',
        #     msd_file=msd_file,
        #     color=WALL_LONG_COLOR,
        #     label='single wall, $t=64\mathrm{{s}}$',
        #     marker='x'
        # ),
        # theory
        dict(
            file='ld_hydro_dpstokes_0.114_L1280',
            source='D0Sk_theory',
            label='theory, no hydrodynamics'
        ),
    ),
    ax = ax,
    plot_against_k=True,
    # allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
    discrete_colors=True,
)
# ax.set_ylim(0.8, 3)

common.save_fig(fig, 'presentations/figures/brennan_test.png')