import common
import numpy as np
import matplotlib.pyplot as plt

# PLOT_AGAINST_K = True
TWO_PI = True
SHOW_TWIN_K_AXIS = False

for file in common.files_from_argv('visualisation/data', 'Ds_from_DDM_'):
    x = 0

    source_names = {
        'DDM': 'DDM',
        'DDM_short': 'DDM short',
        'DDM_long': 'DDM long',
        'f': '$f(k, \Delta t)$',
        'Fs': '$F_s(k, \Delta t)$',
        'f_short': '$f(k, \Delta t)$ short',
        'Fs_short': '$F_s(k, \mathrm{short})$',
        'f_long': '$f(k, \Delta t)$ long',
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
        'DDM': '*',
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

    # timescaleintegral    = 1    * np.array([np.nan,     3.7504e-01, 1.4125e+00, 4.5840e+00, 1.0473e+01, 1.3479e+01, 1.9825e+01, 3.0339e+01, 4.9982e+01, 9.1947e+01, 1.7270e+02, 3.4680e+02, 7.5721e+02, 1.0308e+03, 1.4394e+03, 1.8460e+03])
    # timescale_integral_L = 0.29 * np.array([5.0000e-01, 1.0000e+00, 2.0000e+00, 5.0000e+00, 1.0000e+01, 1.3800e+01, 2.0000e+01, 2.7700e+01, 4.0000e+01, 5.5400e+01, 8.0000e+01, 1.1070e+02, 1.6000e+02, 1.9000e+02, 2.2000e+02, 2.3500e+02])

    # timescaleintegral_x = np.array([	0.103942652329749,0.207885304659498,0.519713261648746,1.03942652329749,1.43440860215054,2.07885304659498,2.87921146953405,4.15770609318996,5.75842293906810,8.31541218637993,11.5064516129032,16.6308243727599,19.7491039426523,22.8673835125448,24.4265232974910])
    # timescaleintegral_y = np.array([0.715084639456640,0.759462918744972,1.46261640042440,2.56073086204352,3.78909538953524,5.41105358248309,6.78255486238759,8.58502158959043,8.95193674486290,9.93854195925673,9.47657037534071,9.06689661448547,9.39221079633647,9.01776477004307,8.02304202449399])

    # for sources in [['f_short'], [f'f_short', 'DDM_short'], ['timescaleint']]:
    # for sources in [['f_short']]:
    
    # this is the stuff that needs to go into cmd31:
    """
    (['boxcounting_collective'             ], False, True, True, True),
    (['boxcounting_collective', 'MSD_short'], False, True, True, True),
    (['f_short', 'MSD_short'], True, True, True, True),
    (['f_short'             ], True, True, True, True),
    (['boxcounting_collective', 'f_short'], False, True, True, False),
    (['boxcounting_collective', 'f_short'], False, False, True, False),
    (['boxcounting_shorttime', 'MSD_short'], False, True, False, False),
    (['boxcounting_shorttime'             ], False, True, False, False),
    """
    
    for (sources, PLOT_AGAINST_K, TWO_PI, logarithmic_y, fix_axes) in [
        (['MSD_short', 'f_short', 'f', 'boxcounting_shorttime'], True, True, True, False),
    ]:
    # for sources in [[f'f_short', 'timescaleint', 'MSD_short', 'MSD_long', 'boxcounting', 'boxcounting_shorttime']]:
        print()

        all_Ds = []
        
        used_sources = []

        figsize = (3.2, 2.8)
        if PLOT_AGAINST_K and SHOW_TWIN_K_AXIS:
            figsize = (3.2, 3.2)
            
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        

        # for source in ['f', 'Fs', 'DDM', 'boxcounting', 'boxcounting_shorttime', 'MSD']:
        # for source in ['boxcounting', 'MSD', 'Fs', 'f', 'DDM']:
        # for source in ['boxcounting', f'Fs_{timescale}', f'f_{timescale}', f'DDM_{timescale}', 'MSD']:
        # for source in ['boxcounting', f'MSD_{timescale}', f'Fs_{timescale}', f'f_{timescale}', f'DDM_{timescale}', f'f_long', f'DDM_long', 'f', 'DDM', 'timescaleint']:
        # for source in [f'f_{timescale}', f'DDM_{timescale}', f'f_long', f'DDM_long', 'f', 'DDM', 'timescaleint']:
        # for source in [f'f_{timescale}', f'DDM_{timescale}', f'DDM_long', 'DDM', 'timescaleint']:
        for source in sources:
            try:
                data = common.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
            except FileNotFoundError:
                print('FileNotFound', source)
                continue
            Ds     = data['Ds']
            D_uncs = data['D_uncs']

            assert np.all(Ds > 0)
            assert np.all(D_uncs > 0)

            if Ds.size == 0:
                print(f'skipping {source} {file}, no Ds found')
                continue

            # pixel_size = None
            source_label = f'{source_names[source]}' if not source.startswith('MSD') else None
            

            if source.startswith('f') or source.startswith('Fs') or source.startswith('DDM'):
                # xs = 1 / data['ks']
                print('max k', data['ks'].max(), '2pi over over D', 2*np.pi/data['ks'].max()/2.82)

                if PLOT_AGAINST_K:
                    xs = data['ks']
                else:
                    if TWO_PI:
                        xs = 2 * np.pi / data['ks']
                        source_label += ' $2\pi/k$'
                    else:
                        xs = 1 / data['ks']
                        source_label += ' $1/k$'
                # print(xs)
                # print(Ds)

                # pixel_size = data['pixel_size']
                # print('HACK')
                # pixel_size = 0.288
                if file == 'eleanorlong':
                    print('REMOVING POINTS!')
                    start = 9
                    xs = xs[start:]
                    Ds = Ds[start:]
                    D_uncs = D_uncs[start:]
                else:
                    print('not')


                # xs = data['ks'] / (2 * np.pi / diameter)
            elif source.startswith('boxcounting'):
                xs = data['Ls']
                if PLOT_AGAINST_K:
                    raise

                # print('THIS IS SO NAUGHTY REMOVE ME')
                # thresh = 0.576
                # Ds = Ds[xs > thresh]
                # D_uncs = D_uncs[xs > thresh]
                # xs = xs[xs > thresh]

            elif source.startswith('timescaleint'):
                if PLOT_AGAINST_K:
                    print('skipping timescaleint')
                    continue
                else:
                    xs = data['Ls']
                    print('tsix', xs, Ds)
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
                xs = np.array([])
                Ds = np.array([])
                D_uncs = np.array([])
            else:
                raise Exception(f'you need to specify the x scale for {source}')


            diameter = data.get('particle_diameter')
            print('must change above line!')
            ten = 10
            if PLOT_AGAINST_K:
                if not np.isnan(diameter):
                    xs *= diameter
                # ten *= diameter
                # if pixel_size:
                #     pixel_size *= diameter
            else:
                if not np.isnan(diameter):
                    xs /= diameter
                # ten /= diameter
                # if pixel_size:
                #     pixel_size /= diameter

            print(xs)

            assert not np.any(np.isnan(xs)), f'nan found in xs from {source}'

            # if not source.startswith('MSD'):
            color = colors.get(source)
            lines = ax.plot(xs, Ds, linestyle='none', marker=markers[source], markersize=4, color=color, label=source_label)
            print(color)
            ax.errorbar(xs, Ds, yerr=D_uncs, linestyle='none', marker='none', alpha=0.6, color=color)

            print('plotting', xs, Ds)

            # if pixel_size:
            #     ax.vlines(pixel_size, Ds.min(), Ds.max(), linestyle='dotted', color='black', label='pixel')
            #     ax.vlines(ten, Ds.min(), Ds.max(), linestyle='dotted', color='grey', label='pixel')

            # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'
            [all_Ds.append(D) for D in Ds]

            used_sources.append(source)


        assert len(all_Ds) > 0, 'no Ds were found at all'
        print('all Ds', all_Ds)

        phi = 0.34
        D_th_self = ( 1 - 1.73*phi )
        D_th_coll = ( 1 + 1.45*phi )
        # ax.hlines([D_th_coll, D_th_self], xs.min(), xs.max(), color='black')

        if logarithmic_y:
            ax.semilogy()

        ylim_expand = 1.2
        if np.nanmax(all_Ds) - np.nanmax(all_Ds) < 0.4:
            ylim_expand = 1.5
        ax.set_ylim(np.nanmin(all_Ds)/ylim_expand, np.nanmax(all_Ds)*ylim_expand)
        if fix_axes:
            ax.set_ylim(0.018, 1)
        ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
        # ax.set_ylabel('$D/D_0$')
        # ax.set_xticks([])
        ax.semilogx()
        if PLOT_AGAINST_K:
            
            if not np.isnan(diameter):
                ax.set_xlabel(r'$k \sigma$')
            else:
                ax.set_xlabel(r'$k$')
            # ax.set_xlim(0.3, 20)

            if SHOW_TWIN_K_AXIS:
                realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
                # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
                realspace_ax.set_xlabel(r'$2\pi/k / \sigma$')
        else:
            ax.set_xlabel(r'$L / \sigma$')
        # ax.set_xlabel(r'$k / (2\pi / \sigma)$')
        # ax.semilogy()
        # ax.legend(fontsize=5, loc='upper left')


        # ax.set_ylim(0.3, 2)
        # print('NORTY YLIM')

        # ax.legend(fontsize=7, loc='lower right')
        ax.legend(fontsize=7)
        # ax.set_title(f'{file}, errorbars not yet all correct')
        name = '_'.join(used_sources)
        filename = f'{file}_{name}'
        if not TWO_PI:
            filename += '_oneoverk'
        # common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/Ds_overlapped_{filename}.pdf', hide_metadata=True)
        common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_{filename}.png', dpi=200)