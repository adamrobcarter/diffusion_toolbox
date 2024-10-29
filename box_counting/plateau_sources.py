import common
import matplotlib.pyplot as plt
import box_counting.D_of_L
import visualisation.Ds_overlapped

for file in common.files_from_argv('box_counting/data', 'counted_'):
    # fig, axs = plt.subplots(5, 2, figsize=(2*3, 5*3))
    fig, axs = plt.subplots(2, 5, figsize=(5*3, 2*3))

    for method_index, method in enumerate(['obs', 'var', 'varmod', 'nmsdfit', 'nmsdfitinter']):
        box_counting.D_of_L.go(file, axs[0][method_index], plateau_source=method, legend_fontsize=5, title=f'plateau:{method}')
        visualisation.Ds_overlapped.go(file, [f'timescaleint_nofit_cropped_{method}', f'timescaleint_{method}', 'MSD_first'], ax=axs[1][method_index], ylim=(0.6, 100))

    common.save_fig(fig, f'box_counting/figures_png/plateau_sources_{file}.png', dpi=200)