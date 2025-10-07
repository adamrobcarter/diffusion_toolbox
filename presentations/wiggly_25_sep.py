import common
import visualisation.Ds_overlapped_mult
import isf.show_S_of_k
import matplotlib.pyplot as plt
import numpy as np
import isf.show_Fs_overlayed
import matplotlib.image
import matplotlib.offsetbox

figs = '/home/acarter/presentations/wiggly_25_sep/figures'

WALL_LONG_COLOR = 'royalblue'
WALL_SHORT_COLOR = 'cornflowerblue'
OPEN_LONG_COLOR = 'indianred'
OPEN_SHORT_COLOR = 'lightcoral'

def show_png_on_axis(ax, file, coords, size, color='black'):
    image = matplotlib.image.imread(file).copy()
    print(file, image.dtype, image.min(), image.max(), image.shape)
    imagebox = matplotlib.offsetbox.OffsetImage(image, zoom=size)   
    imagebox.set_data(image)

    ab = matplotlib.offsetbox.AnnotationBbox(
        imagebox,
        coords,        # position in DATA units
        xycoords="data", # <-- use data coordinates
        # frameon=False,
    )
    
    frame = ab.patch
    frame.set_edgecolor(color)     # frame/border color
    # frame.set_linewidth(2)         # thickness

    ax.add_artist(ab)

    return imagebox, ab



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
    show_twin_k_axis=True,
    allow_rescale_x=True
)

common.save_fig(fig, f'{figs}/collective_diffusion_1.pdf', hide_metadata=True)


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
    show_twin_k_axis=True,
    allow_rescale_x=True
)
common.add_exponential_index_indicator(ax, -1, (1e-1, 1e1), 'k', x_limits=(1e-2, 2e-1), )
ib1, ab1 = show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (1e1, 5), 0.04, color=OPEN_LONG_COLOR)


common.save_fig(fig, f'{figs}/collective_diffusion_2.pdf', hide_metadata=True)



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
    show_twin_k_axis=True,
    allow_rescale_x=True
)
ib2, ab2 = show_png_on_axis(ax, 'presentations/libmobility/side_view_2.png', (1e1, 2.5), 0.04, color=WALL_LONG_COLOR)
# show_png_on_axis(ax, 'presentations/libmobility/side_view_potential.png', (1e1, 2.5), 0.04, color=WALL_LONG_COLOR)


common.save_fig(fig, f'{figs}/collective_diffusion_3.pdf', hide_metadata=True)