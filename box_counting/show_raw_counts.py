import common
import matplotlib.pyplot as plt
import matplotlib.cm

for file in common.files_from_argv('box_counting/data', 'raw_counts_'):
    data = common.load(f'box_counting/data/raw_counts_{file}.npz')
    counts = data['counts']
    counts = counts[0] # box_size_index == 0

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.2))

    MAX_TO_PLOT = 10
    num_plotted = 0

    for x in range(len(counts)):
        for y in range(len(counts[x])):
            # ax.plot(counts[x][y])
            chosen = num_plotted == 3
            color = matplotlib.cm.afmhot(0.3) if chosen else matplotlib.cm.Greys(num_plotted/MAX_TO_PLOT)
            zorder = 2 if chosen else 0
            linewidth = 3 if chosen else 1
            ax.stairs(counts[x][y], color=color, zorder=zorder, linewidth=linewidth, alpha=0.7)

            num_plotted += 1
            if num_plotted > MAX_TO_PLOT:
                break

    ax.set_xlim(0, 200)
    ax.set_xlabel('$t$ (s)')
    ax.set_ylabel('$N(t)$')

    # common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/raw_counts_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'box_counting/figures_png/raw_counts_{file}.png')