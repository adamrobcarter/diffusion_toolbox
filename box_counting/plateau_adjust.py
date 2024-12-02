import matplotlib.pyplot as plt
import common
import box_counting.D_of_L
import visualisation.Ds_overlapped_mult

SHOW_THEORY = True
SHOW_TIMESCALEINTEGRAL_FIT = True
DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE = True
SHOW_LEGEND = True

LATE_CN_ALPHA = 0.1
LABELS_ON_PLOT = False
LABELS_ON_PLOT_Y_SHIFT = 1.2
LABELS_ON_PLOT_X_SHIFT = 1.3

MAX_NUM_BOXES = 10

UNC_INCLUDES_NMSD_UNC = False

NOFIT_CROP_THRESHOLD = 1e-4

RESCALE_X_L2 = 1
# RESCALE_X = RESCALE_X_L2
RESCALE_X = None

PLATEAU_SOURCE = 'var'


if __name__ == '__main__':  
    for file in common.files_from_argv('box_counting/data', 'counted_'):

        adjusts = 0.98, 0.99, 1.0, 1.01, 1.02
        # adjusts = -0.01, -0.001, 0, 0.001, 0.01
        for adjust in adjusts:
            dead_fig, dead_axs = plt.subplots(1, 1, figsize=(6, 5))

            filelabel = str(adjust).replace('-', 'm')
            box_counting.D_of_L.go(file, ax=dead_axs, plateau_source=PLATEAU_SOURCE, show_theory=SHOW_THEORY,
                labels_on_plot=LABELS_ON_PLOT, max_num_boxes=MAX_NUM_BOXES, show_short_fits=SHOW_TIMESCALEINTEGRAL_FIT, show_long_fits=SHOW_TIMESCALEINTEGRAL_FIT,
                save_data=True, show_legend=SHOW_LEGEND, rescale_x=RESCALE_X,
                plateau_adjust=adjust,
                # plateau_offset=adjust,
                plateau_source_suffix=f'_pa{filelabel}',
            )

        fig, ax = plt.subplots(1, 1)
        visualisation.Ds_overlapped_mult.go(
            files=[file],
            ax=ax,
            sources=[f'timescaleint_fixexponent_{PLATEAU_SOURCE}_pa' + str(a).replace('-', 'm') for a in adjusts] + ['D_of_L_theory'],
            colors=[[common.colormap(i, 0, len(adjusts)) for i in range(len(adjusts))] + ['black']],
            legend_fontsize=7,
            linestyles='-',
        )
        ax.set_ylim(0.5, 3)
        common.save_fig(fig, f'box_counting/data/Ds_overlapped_plateau_adjust_{file}.png')