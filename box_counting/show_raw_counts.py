import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

MAX_TIME = 200
MAX_TO_PLOT = 1
LABELS_ON_PLOT_Y_ADD = 0.4

def go(file, ax):
    data = common.load(f'box_counting/data/counted_counts_{file}.npz')
    all_counts = data['counts']
    box_sizes = data['box_sizes']
    sigma = data.get('particle_diameter')

    box_size_indices = [6, 13] # 6, 13 is nice
    box_size_indices = [7, 11]

    box = [0, 5]

    i = -1

    for box_size_index in box_size_indices:
        counts = all_counts[box_size_index]
        print(f'doing L={box_sizes[box_size_index]:.1f}')
        print('avg', counts[:, :, :MAX_TIME].mean())
        
        color = common.colormap(box_size_index, 0, len(box_sizes))
        num_plotted = 0
        break_outer = False

        i += 1

        for x in range(box[i], counts.shape[0]):
            for y in range(counts.shape[1]):
                # ax.plot(counts[x][y])
                chosen = num_plotted == 0
                # color = common.colormap(0.5) if chosen else matplotlib.cm.Greys(np.interp(num_plotted/MAX_TO_PLOT, (0, 1), (0, 0.5)))
                zorder = 2 if chosen else 0
                linewidth = 2 if chosen else 1
                # print(counts[x, y, :])
                ax.stairs(counts[x, y, :MAX_TIME]+0.02, baseline=counts[x, y, 0], color=color, zorder=zorder, linewidth=linewidth, alpha=1)
                # +0.1 just visually helps it pop a bit at N(t) == 0
                print('thi', counts[x, y, :MAX_TIME].mean())

                
                t_index_for_text = 0

                if sigma and not np.isnan(sigma):
                    L_label = rf'$L={box_sizes[box_size_index]/sigma:.2g}\sigma$'
                else:
                    L_label = rf'$L={box_sizes[box_size_index]:.2g}$'

                if counts[x, y, :50].max() == 0 or counts[x, y, :50].min() == 0:
                    y = counts[x, y, :50].max() + LABELS_ON_PLOT_Y_ADD
                    va = 'baseline'
                else:
                    y = counts[x, y, :50].min() - LABELS_ON_PLOT_Y_ADD
                    va = 'top'
                
                ax.text(5, y, L_label,
                        ha='left', va=va, color=color, fontsize=12)
            

                num_plotted += 1
                if num_plotted >= MAX_TO_PLOT:
                    break_outer = True
                    break
            if break_outer:
                break
            

    ax.set_xlim(0, MAX_TIME)
    ax.set_xlabel('$t$ (s)')
    ax.set_ylabel('$N(t)$')


if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data', 'counted_counts_'):
        fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.2))

        go(file, ax)
        # common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/raw_counts_{file}.pdf', hide_metadata=True)
        common.save_fig(fig, f'box_counting/figures_png/raw_counts_{file}.png')