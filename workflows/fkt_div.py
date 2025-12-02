
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
import numpy as np
import matplotlib.pyplot as plt
import visualisation.Ds_overlapped_mult
import common

fig, ax = plt.subplots(figsize=(5,5))

L = np.array([361.6])
get_L = lambda c : rf'$L_x={L[0]*c:.0f}\mathrm{{\mu m}}$'
visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file = 'eleanorlong001_crop1.0',
            source = 'f_first_first',
        ),
        dict(
            file = 'eleanorlong001_crop0.5',
            source = 'f_first_first',
        ),
        dict(
            file = 'eleanorlong001_crop0.25',
            source = 'f_first_first',
        ),
        dict(
            file = 'eleanorlong001_crop0.125',
            source = 'f_first_first',
        )
    ],
    ax=ax,
    plot_against_k=True,
    show_twin_k_axis=False,
    # legend_fontsize=LEGEND_FONTSIZE,
    # markers='d',
    # file_labels=[get_L(1), get_L(0.5), get_L(0.25), get_L(0.125)],
    # source_labels=['', ''],
    # colors=[[COLORMAP_FKT_CROPS(0.0)], [COLORMAP_FKT_CROPS(0.33)], [COLORMAP_FKT_CROPS(0.66)], [COLORMAP_FKT_CROPS(1.0)]]
)
ax.semilogy()
common.add_exponential_index_indicator(ax, exponent=-2, anchor=(0.03, 10), xlabel='k')
common.save_fig(fig, 'workflows/figures/fkt_div.png')