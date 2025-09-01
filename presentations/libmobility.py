import common
import visualisation.Ds_overlapped_mult
import isf.show_S_of_k
import matplotlib.pyplot as plt
import numpy as np
import isf.show_Fs_overlayed

figs = 'presentations/libmobility'

fig, axs = plt.subplot_mosaic(
    """
    ac
    bc
    """,
    gridspec_kw=dict(width_ratios=(1, 2)),
    figsize=(6, 4)
)

visualisation.Ds_overlapped_mult.go(
    (
        ('ld_hydro_nbody_0.114_L2560_t450_1', 'f_t1'),
        ('ld_hydro_nbody_0.114_L2560', 'f_t1024'),
        ('ld_hydro_nbody_0.114_L2560_t450_1', 'D0Sk_theory'),
        ('ld_hydro_nbody_open_0.114_L2560', 'f_t1024'),
    ),
    ax = axs['c'],
    plot_against_k=True,
    allow_rescale_y=False,
    # labels=('hydro above')
    legend_fontsize=7,
    show_twin_k_axis=False,
)
# ax.set_ylim(0.5, 4)
axs['c'].set_ylim(3e-2, 2e0)
axs['c'].set_xlim(2.1e-2, 1e1)
axs['c'].semilogy()
common.add_exponential_index_indicator(axs['c'], -1, (1e-1, 6e-1), 'k', x_limits=(1e-2, 3e-1))
# axs['c'].grid()

isf.show_S_of_k.go(
    'ld_hydro_nbody_0.114_L2560',
    axs['a'],
    source='F',
    show_realspace_axis=False
)
axs['a'].set_ylim(0.55, 1.15)


isf.show_Fs_overlayed.go(
    'ld_hydro_nbody_0.114_L2560',
    file_i=0,
    ax=axs['b'],
    target_ks=(10, 1, 0.1, 0.01)
)



common.save_fig(fig, f'{figs}/monolayer.png')