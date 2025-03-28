import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common

path = '/home/acarter/presentations/wiggly_25_march/figures'

for i in range(2, 6):
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    visualisation.Ds_overlapped_mult.go(
        colors=[
            ['black',
            common.colormap(0),
            common.colormap(0.3),
            common.colormap(0.6),
            common.colormap(0.9)][:i]
        ],
        ax      = ax,
        files   = ['sim_nohydro_011_L640'], 
        sources = [
            'D0Sk_theory',
            'f_t0.5',
            'f_t2',
            'f_t8',
            'f_t32',
        ][:i],
        legend_fontsize=8
    )
    ax.set_ylim(0.8, 2.5)
    
    common.save_fig(fig, f'{path}/mixt{i}.pdf', hide_metadata=True)

DS_HYDRO_YLIM = (0.7, 4)

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
visualisation.Ds_overlapped_mult.go(
    ax      = ax,
    files   = ['sim_nohydro_011_L640_longer_mergedD', 'sim_hydro_011_L640_longer_mergedD'], 
    sources = [
        'D0Sk_theory',
        'f_first_first',
    ],
    file_labels = ['nohydro', 'hydro'],
    legend_fontsize=8
)
ax.set_ylim(*DS_HYDRO_YLIM)
common.save_fig(fig, f'{path}/hydro_fkt.pdf', hide_metadata=True)

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
visualisation.Ds_overlapped_mult.go(
    ax      = ax,
    files   = ['sim_nohydro_011_L640_longer_mergedD', 'sim_hydro_011_L640_longer_mergedD'], 
    sources = [
        'D_of_L_theory',
        'timescaleint_nofit_cropped_var',
    ],
    file_labels = ['nohydro', 'hydro'],
    legend_fontsize=8
)
ax.set_ylim(*DS_HYDRO_YLIM)
common.save_fig(fig, f'{path}/hydro_counting.pdf', hide_metadata=True)