import common
import matplotlib.pyplot as plt
import box_counting.D_of_L
import box_counting.msd_single
import visualisation.Ds_overlapped

sources = ['var', 'varmod', 'nmsdfit', 'nmsdfitinter', 'sDFT']

for file in common.files_from_argv('box_counting/data', 'counted_'):
    # fig, axs = plt.subplots(5, 2, figsize=(2*3, 5*3))
    fig, axs = plt.subplots(2, len(sources), figsize=(len(sources)*3, 2*3))

    dead_fig, dead_ax = plt.subplots(1, 1)

    for method_index, method in enumerate(sources):
        print()
        print('doing msd_single')
        box_counting.msd_single.go(
            file,
            ax=dead_ax,
            timescaleint_replacement_plateau_source=method,
        )
        print()
        print('doing D_of_L')
        box_counting.D_of_L.go(
            file,
            ax=axs[0][method_index],
            plateau_source=method,
            legend_fontsize=5,
            title=f'plateau:{method}',
            save_data=True,
        )
        print()
        print('doing Ds_overlapped')
        visualisation.Ds_overlapped.go(
            file,
            [f'boxcounting_collective_{method}', f'timescaleint_nofit_cropped_{method}', f'timescaleint_{method}', 'MSD_first', 'D_of_L_theory'],
            ax=axs[1][method_index],
            ylim=(0.2, 50),
            legend_fontsize=5,
        )

    common.save_fig(fig, f'box_counting/figures_png/plateau_sources_{file}.png', dpi=200)