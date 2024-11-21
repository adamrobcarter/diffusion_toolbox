import common
import matplotlib.pyplot as plt
import box_counting.D_of_L
import visualisation.Ds_overlapped

sources = ['var', 'varmod', 'nmsdfit', 'nmsdfitinter']

for file in common.files_from_argv('box_counting/data', 'counted_'):
    # fig, axs = plt.subplots(5, 2, figsize=(2*3, 5*3))
    fig, axs = plt.subplots(2, len(sources), figsize=(len(sources)*3, 2*3))

    for method_index, method in enumerate(sources):
        box_counting.D_of_L.go(
            file,
            ax=axs[0][method_index],
            plateau_source=method,
            legend_fontsize=5,
            title=f'plateau:{method}',
        )
        visualisation.Ds_overlapped.go(
            file,
            [f'timescaleint_nofit_cropped_{method}', f'timescaleint_{method}', 'MSD_first', 'D_of_L_theory'],
            ax=axs[1][method_index],
            ylim=(0.2, 50),
            legend_fontsize=5,
        )

    common.save_fig(fig, f'box_counting/figures_png/plateau_sources_{file}.png', dpi=200)