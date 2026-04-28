import matplotlib.pyplot as plt
import workflows.thesis.common
import common
import visualisation.Ds_overlapped_mult
import numpy as np

fig, ax = plt.subplots(figsize=workflows.thesis.common.figsize_small,)


SIMULATION_SOURCE_011 = 'sim_nohydro_011_L1280'
SIMULATION_SOURCE_002 = 'sim_nohydro_002_L1280'
sources = [
    'timescaleint_nofit_cropped_var',
    'boxcounting_collective_var',
    'f_first_first',
]
MERGE_TYPE_NMSD_FIT = 'merged' # or {MERGE_TYPE}
MERGE_TYPE_TIMESCALEINT = 'mergedD' # or {MERGE_TYPE}

fig8_props = dict(
    logarithmic_y=False,
    # legend_fontsize=LEGEND_FONTSIZE,
    file_labels=[''],
    # sources=      ['timescaleint_nmsdfitinter', 'D_of_L_theory',      'f_first_first', 'D0Sk_theory'],
    sources=      ['timescaleint_nofit_cropped_var',          'D_of_L_theory',      'f_first_first', 'D0Sk_theory'],
    # colors=       [[COLOR_PHI011,               COLOR_PHI011_THEORY,  COLOR_PHI011_alt,    COLOR_PHI011_THEORY]],
    # markers=      [[MARKER_COUNTING,            'none',               MARKER_FKT,      'none', ]],
    # linestyles=   [['none',                     LINESTYLE_COUNTING,   'none',          LINESTYLE_FKT]],
    source_labels=['Countoscope',               'Countoscope theory', '$f(k, t)$',     '$f(k, t)$ theory']
)

visualisation.Ds_overlapped_mult.go(
    [
        dict(
            file = 'sim_nohydro_011_L1280_longer_mergedD',
            source = 'f_first_first',
            plot_index = np.index_exp[10:],
            marker = 'o',
        ),
        dict(
            file = 'sim_nohydro_011_L1280_longer_mergedD',
            source = 'D0Sk_theory',
            sigma = 2.97,
            phi = 0.114,
            msd = 'sim_nohydro_011_L1280_longer_mergedD'
        )
    ],
    ax=ax,
    allow_rescale_x = True
    # **fig8_props
)
ax.get_legend().remove()
# ax_a.set_ylim(*DS_OVERLAPPED_YLIM)
# ax_a.yaxis.set_major_locator(ticks_0p5)
ax.set_xlim(1.3e-1, 3.5e1)


common.save_fig(fig, f'workflows/thesis/figures/intro_Dk.pdf', hide_metadata=True)