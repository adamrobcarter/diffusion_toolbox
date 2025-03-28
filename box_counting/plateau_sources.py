import common
import matplotlib.pyplot as plt
import box_counting.D_of_L
import box_counting.msd_single
import visualisation.Ds_overlapped_mult
from box_counting.msd_single import PLATEAU_SOURCES, PLATEAU_SOURCE_NAMES

def go(file, sources, axs, Ds_overlapped_kwargs={}, D_of_L_kwargs={}, file_for_D=None, file_CN_lowL=None, file_CN_highL=None):
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

        all_D_of_L_kwargs = dict(
            ax=axs[0][method_index],
            plateau_source=method,
            legend_fontsize=5,
            save_data=True,
            disable_ylabel=False if method_index == 0 else True,
            **D_of_L_kwargs,
        )

        if True:
            box_counting.D_of_L.go(
                f'{file}',
                box_size_indices=[0, 3, 6, 9, 12],
                **all_D_of_L_kwargs
            )
            box_counting.D_of_L.go(
                f'{file}_longer',
                box_size_indices=[15, 18, 21, 24, 27],
                **all_D_of_L_kwargs
            )
        else:
            box_counting.D_of_L.go(
                f'{file}_longer_merged',
                **all_D_of_L_kwargs
            )
        axs[0][method_index].set_title(f'Plateau: {PLATEAU_SOURCE_NAMES[method]}', fontsize=10)
        axs[0][method_index].set_xlim(0.5, 0.8e5)
        print()
        print('doing Ds_overlapped')
        visualisation.Ds_overlapped_mult.go(
            files   = [f'{file}_longer_merged', f'{file}_longer_mergedD'],
            sources = [
                f'boxcounting_collective_{method}',
                f'timescaleint_nofit_cropped_{method}',
                'D_of_L_theory'
            ],
            ax=axs[1][method_index],
            legend_fontsize=7,
            file_labels = ['', ''],
            disable_ylabel=False if method_index == 0 else True,
            **Ds_overlapped_kwargs
        )
        axs[1][method_index].set_yscale('linear')
        # if file == 'eleanorlong001':
        #     axs[1][method_index].set_ylim(0.9*0.04, 1*0.07)
        # else:
        #     axs[1][method_index].set_ylim(0.9*0.04, 3*0.04)
        # break

if __name__ == '__main__':
    
    for file in common.files_from_argv('box_counting/data', 'counted_'):
        
        # sources = [
        #     'var',
        #     'histogram'
        # ]
        sources = ['var', 'varmod', 'nmsdfit', 'sDFT']
        fig, axs = plt.subplots(2, len(sources), figsize=(len(sources)*3, 2*3))

        go(file,
           sources=sources,
           axs=axs,
           
        )

        common.save_fig(fig, f'box_counting/figures_png/plateau_sources_{file}.png', dpi=200)