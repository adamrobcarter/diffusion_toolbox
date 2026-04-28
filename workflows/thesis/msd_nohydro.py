import MSD.show_D
import MSD.show
import workflows.thesis.common
import matplotlib.pyplot as plt
import common
import numpy as np

files = [
    'ld_self_flat_0.2_singlewall_L1280_t2m_20ms_dt20_nolub_unwrap_2d',
    'ld_self_flat_0.2_singlewall_L1280_t42m_1s_dt20_nolub_unwrap_2d',
    # 'ld_self_flat_0.1_singlewall_L1280_t2m_20ms_dt20_nolub_unwrap_2d',
]
index_exps = [
    np.index_exp[2:-4400],
    np.index_exp[10:-1000],
    np.index_exp[:]
]

fig, ax = plt.subplots(1, 1, figsize=workflows.thesis.common.figsize_halfpage)

for i in range(len(files)):
    MSD.show.go(
        files[i],
        ax=ax,
        color=workflows.thesis.common.discrete_colors[0],
        index_exp=index_exps[i],
        markersize=5,
        SHOW_SHORT_FIT=False,
        SHOW_LONG_FIT=False,
    )
ax.autoscale()

common.save_fig(fig, 'workflows/thesis/figures/msd_nohydro.pdf', hide_metadata=True)

fig, ax = plt.subplots(1, 1, figsize=workflows.thesis.common.figsize_halfpage)

for i in range(len(files)):
    MSD.show_D.go(
        files[i],
        ax=ax,
        color=workflows.thesis.common.discrete_colors[0],
        index_exp=index_exps[i],
        markersize=5,
    )

# MSD.show_D.go('ld_self_flat_0.2_singlewall_L1280_t1s_0ms_dt0_nolub_unwrap_2d', ax=ax)
# MSD.show_D.go('ld_self_flat_0.2_singlewall_L1280_t42m_1s_dt10_nolub_unwrap_2d', ax=ax)

common.save_fig(fig, 'workflows/thesis/figures/D_from_msd_nohydro.pdf', hide_metadata=True)