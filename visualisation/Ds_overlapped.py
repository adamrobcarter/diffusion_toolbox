import common
import numpy as np
import matplotlib.pyplot as plt

for file in common.files_from_argv('visualisation/data', 'Ds_from_DDM_'):
    x = 0

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

    for timescale in ['short', 'long']:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

        # for source in ['f', 'Fs', 'DDM', 'boxcounting', 'boxcounting_shorttime', 'MSD']:
        # for source in ['boxcounting', 'MSD', 'Fs', 'f', 'DDM']:
        for source in ['boxcounting', f'Fs_{timescale}', f'f_{timescale}', f'DDM_{timescale}', 'MSD']:
            data = common.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
            Ds     = data['Ds']
            D_uncs = data['D_uncs']

            if source.startswith('f') or source.startswith('Fs') or source.startswith('DDM'):
                xs = 2 * np.pi / data['ks']
            elif source.startswith('boxcounting'):
                xs = data['Ls']
            elif source == 'MSD':
                xs = ax.get_xlim()[0]
            else:
                raise Exception('you need to specify the x scale for this type')
            
            # if source in 

            ax.errorbar(xs, Ds, D_uncs, linestyle='none', marker='_', label=f'$D$ from {source_names[source]}')

            [all_Ds.append(D) for D in Ds]
        ax.set_ylim(-0.12, 0.3)

        # ax.set_ylim(0, np.median(all_Ds)*2)
        ax.set_ylabel('$D$')
        ax.set_xticks([])
        ax.semilogx()
        ax.set_xlabel('length scale ($\mathrm{\mu m}$)')
        # ax.semilogy()
        # fig.legend()
        fig.legend(loc='upper right', bbox_to_anchor=(0.96, 0.4), fontsize=8)
        # ax.set_title(f'{file}, errorbars not yet all correct')
        common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_{timescale}_{file}.png', dpi=200)