import common
import numpy as np
import matplotlib.pyplot as plt
import sys
import DDM.show
import matplotlib.cm

# INDEX = 5

source_names = {
    'DDM': 'DDM',
    'DDM_short': 'DDM short',
    'DDM_long': 'DDM long',
    'f': '$f(k, t)$',
    'Fs': '$F_s(k, t)$',
    'f_short': '$f(k, t)$',
    'Fs_short': '$F_s(k, \mathrm{short})$',
    'f_long': '$f(k, t)$ long',
    'Fs_long': '$F_s(k, \mathrm{long})$',
    'boxcounting': 'counting full fit',
    'MSD': 'MSD',
    'MSD_short': 'MSD',
    'MSD_long': 'MSD long time',
    'MSD_first': 'MSD first',
    'boxcounting_shorttime': 'Countoscope short time fit',
    'boxcounting_first_quad': 'Countoscope short time',
    'boxcounting_collective': 'Countoscope full fit',
    'timescaleint': 'timescale integral'
}

colors = {
    'DDM_short': 'tab:purple',
    'f_short': 'tab:green',
    'f_long': 'tab:green',
    'Fs_short': 'tab:green',
    # 'boxcounting': 'counting',
    'MSD_short': 'tab:blue',
    'timescaleint': 'tab:green',
    'boxcounting_collective': 'tab:orange',
    'boxcounting_shorttime': 'tab:orange',
    'boxcounting_first_quad': 'tab:orange',
}

markers = {
    'DDM_short': '*',
    'f_short': 'x',
    'DDM_long': '*',
    'f_long': 'x',
    'DDM': '_',
    'f': 'x',
    'MSD': '_',
    'MSD_short': '_',
    'MSD_long': '_',
    'MSD_first': '_',
    'boxcounting': '_',
    'boxcounting_first_quad': '_',
    'boxcounting_shorttime': '_',
    'boxcounting_collective': '_',
    'timescaleint': 'o'
}

for (sources,) in [
    (['DDM'],),
    # (['f'],),
    # (['f_short'],),
]:
        
    figsize = (5.2, 4.8)
        
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    x = 0

    files = []
    xtick_pos = []
        
    for file in common.files_from_argv('visualisation/data', 'Ds_from_DDM_'):


        x += 1
        print()

        all_Ds = []

        for source in sources:

            try:
                data = common.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
            except FileNotFoundError:
                print('FileNotFound', source)
                continue


            Ds     = data['Ds']
            D_uncs = data['D_uncs']
            NAME   = str(data.get('NAME', file))

            if Ds.size == 0:
                print(f'skipping {source} {file}, no Ds found')
                continue

            # pixel_size = None
            source_label = f'{source_names[source]}' if not source.startswith('MSD') else None
            

            if source.startswith('f') or source.startswith('Fs') or source.startswith('DDM'):
                
                # xs = 1 / data['ks']
                print('max k', data['ks'].max(), '2pi over over D', 2*np.pi/data['ks'].max()/2.82)

                color = data['ks']
                # print(xs)
                # print(Ds)

                # pixel_size = data['pixel_size']
                # print('HACK')
                # pixel_size = 0.288
                # if file == 'eleanorlong':
                #     print('REMOVING POINTS!')
                #     start = 9
                #     color = color[start:]
                #     Ds = Ds[start:]
                #     D_uncs = D_uncs[start:]
                # else:
                #     print('not')


                # xs = data['ks'] / (2 * np.pi / diameter)
            elif source.startswith('boxcounting'):
                color = data['Ls']

            elif source.startswith('timescaleint'):
                color = data['Ls']
                print('tsix', color, Ds)
            elif source.startswith('MSD'):
                # xs = ax.get_xlim()[0]+0.1
                pixel_size = []
                # print(10/(Ds/0.0455840), 'D')
                ax.hlines(Ds, *ax.get_xlim(), color=colors['MSD_short'], linestyle='dotted', label=f'{source_names[source]}')
                # ax.hlines(Ds/0.0455840 / common.structure_factor_2d_hard_spheres(0.001, 0.34, diameter), 0.04, 25, color='grey', linestyle='dotted')
                # ax.hlines(Ds/0.0455840 / common.structure_factor_2d_hard_spheres(0.01, 0.34, diameter), 0.04, 25, color='grey', linestyle='dotted')
                # ax.hlines(Ds/0.0455840 / common.structure_factor_2d_hard_spheres(0.1, 0.34, diameter), 0.04, 25, color='grey', linestyle='dotted')
                # print(common.structure_factor_2d_hard_spheres(0.1, 0.34, diameter))
                # print(common.structure_factor_2d_hard_spheres(0.01, 0.34, diameter))
                # print(common.structure_factor_2d_hard_spheres(0.001, 0.34, diameter))
                # continue
                print('MSD errors hacked')
                ax.fill_between(ax.get_xlim(), Ds[0]*0.97, Ds[0]*1.03, facecolor=colors['MSD_short'], alpha=0.5)
                print('msd unc', D_uncs)
                all_Ds.append(Ds[0])
                color = np.array([])
                Ds = np.array([])
                D_uncs = np.array([])
            else:
                raise Exception(f'you need to specify the color scale for {source}')


            # plot_color = colors.get(source)
            xs = np.full_like(Ds, x)

            # xs = xs[[INDEX]]
            # Ds = Ds[[INDEX]]
            # D_uncs = D_uncs[[INDEX]]

            # make transparency proportional to error
            rel_error = D_uncs/Ds # in (0, DDM.show.D_ERROR_THRES)
            error0to1 = np.interp(rel_error, (0, DDM.show.D_ERROR_THRESH), (1, 0))
            
            if channel := data.get('channel'):
                if channel == 'red':
                    cmap = matplotlib.cm.Reds
                    color2 = 'red'
                elif channel == 'green':
                    cmap = matplotlib.cm.Greens
                    color2 = 'green'
            else:
                cmap = None

            assert not np.any(np.isnan(xs))

            print(cmap(color))

            for i in range(xs.size):
                ax.errorbar(xs[i], Ds[i], yerr=D_uncs[i], color=cmap(np.interp(color[i], (color.min(), color.max()), (0.2, 1))), marker=markers[source], label=source_label)
            
            
            files.append(NAME)
            xtick_pos.append(x)
            # ax.errorbar(color, Ds, yerr=D_uncs, linestyle='none', marker='none', alpha=0.6)

            # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'
            [all_Ds.append(D) for D in Ds]

    # fig.colorbar(scatter)

    
    # plt.colorbar.ColorbarBase(ax, cmap=cmap)#, orientation='vertical')   

    assert len(all_Ds) > 0, 'no Ds were found at all'

    ylim_expand = 1.2
    if np.nanmax(all_Ds) - np.nanmax(all_Ds) < 0.4:
        ylim_expand = 1.5
    # ax.set_ylim(np.nanmin(all_Ds)/ylim_expand, np.nanmax(all_Ds)*ylim_expand)
    ax.set_ylim(0, np.nanquantile(all_Ds, 0.9))


    ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
    # ax.set_ylabel('$D/D_0$')
    ax.set_xticks(xtick_pos, files, rotation=-90)

    # ax.semilogy()
    
    # ax.legend(fontsize=7)
    # ax.set_title(f'{file}, errorbars not yet all correct')

    sources_name = '_'.join(sources)
    filename = 'Ds_mult_' + sources_name + '_' + '_'.join(sys.argv[1:])
    
    common.save_fig(fig, f'visualisation/figures_png/{filename}.png', dpi=200)