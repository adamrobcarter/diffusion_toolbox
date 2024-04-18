import common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

for file in common.files_from_argv('visualisation/data', 'Ds_from_DDM_'):
    x = 0
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))

    source_names = {
        'DDM': 'DDM',
        'f': '$f(k, t)$',
        'Fs': '$F_s(k, t)$',
        'f_short': '$f(k, \mathrm{short})$',
        'Fs_short': '$F_s(k, \mathrm{short})$',
        'f_long': '$f(k, \mathrm{long})$',
        'Fs_long': '$F_s(k, \mathrm{long})$',
        'boxcounting': 'counting',
        'MSD': 'MSD',
        'boxcounting_shorttime': 'Box Counting s.t.a.'
    }

    all_Ds = []

    arrowname = {
        'DDM': '$k$',
        'f': '$k$',
        'Fs_short': '$k$',
        'f_short': '$k$',
        'Fs_long': '$k$',
        'f_long': '$k$',
        'Fs': '$k$',
        'boxcounting': '$L$',
        'MSD': None,
        'boxcounting_shorttime': '$L$'
    }

    # for source in ['f', 'Fs', 'DDM', 'boxcounting', 'boxcounting_shorttime', 'MSD']:
    # for source in ['boxcounting', 'MSD', 'Fs', 'f', 'DDM']:
    for source in ['boxcounting', 'MSD', 'Fs_short', 'Fs_long', 'f_short', 'f_long', 'DDM']:
        data = np.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
        Ds     = data['Ds']
        D_uncs = data['D_uncs']

        xs = x + np.linspace(-0.4, 0.4, num=Ds.size)
        if Ds.size == 1:
            xs = [x]

        label = f'$D$ from {source_names[source]}'
        label = f'{source_names[source]}'
        plotted = ax.errorbar(xs, Ds, D_uncs, linestyle='none', marker='_')

        # label_y = Ds.min()-0.01
        if file == 'eleanor0.01':
            label_y = 0.035
            mult = 1
        if file == 'eleanor0.34':
            pass
        label_y = 0
        mult = 4

        if arrowname[source]:
            # ax.arrow(xs[0], Ds.min()-0.01, xs[-1]-xs[0], 0, color=plotted[0].get_color(), width=0.0001)
            ax.annotate("", xy=(xs[-1], label_y-0.0), xytext=(xs[0], label_y-0.0), arrowprops=dict(arrowstyle="->", color=plotted[0].get_color(), linewidth=0.7))
            ax.text((xs[-1]+xs[0])/2, label_y-0.005*mult, arrowname[source], color=plotted[0].get_color(), fontsize=8, ha='center')

        ax.text((xs[-1]+xs[0])/2, label_y-0.015*mult, label, color=plotted[0].get_color(), fontsize=10, rotation=-45, rotation_mode='anchor', horizontalalignment='center', verticalalignment='center')

        x += 1

        [all_Ds.append(D) for D in Ds]

    ax.relim() # tell mpl to ignore errorbars when
    ax.autoscale_view() # calcing axis limits

    ax.set_ylim(0, ax.get_ylim()[1])

    if file == 'eleanor0.01':
        ax.set_ylim(0, np.median(all_Ds)*2)
    if file == 'eleanor0.34':
        pass
    ax.set_ylim(-0.12, 0.3)
    # ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'{x:0.2f}' if x >= 0 else ''))
    ax.set_yticks(np.arange(0, 0.3, 0.05))
    # ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1]+0.0)
    # ax.hlines(np.median(all_Ds) * (1 - 0.34))
    ax.set_ylabel('$D$ ($\mathrm{\mu m}^{-1}$)')
    ax.set_xticks([])
    # fig.legend()
    # fig.legend(loc='upper right', bbox_to_anchor=(0.96, 0.9), fontsize=9)
    # ax.set_title(f'{file}')

    common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/Ds_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'visualisation/figures_png/Ds_{file}.png', dpi=200)