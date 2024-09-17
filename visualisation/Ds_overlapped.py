import common
import numpy as np
import matplotlib.pyplot as plt
import sys

SHOW_TWIN_K_AXIS = False
PRESENT_SMALL = False

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
    'MSD_centre_of_mass_onepoint': 'MSD CoM',
    'MSD_centre_of_mass_proximity': 'MSD CoM prox',
    'boxcounting_shorttime': 'Countoscope short time fit',
    'boxcounting_first_quad': 'Countoscope short time',
    'boxcounting_collective': 'Countoscope full fit',
    'timescaleint': 'timescale integral',
    'timescaleint_nofit': 'timescale integral (no fit)',
    'D0Sk': r'$D_{MSD}/S(k)$',
    'C_N_simplefit': '$C_N$ fit',
}

colors = {
    'DDM_short': 'tab:purple',
    'f': 'lime',
    'f_short': 'tab:green',
    'f_long': 'yellowgreen',
    'Fs_short': 'tab:green',
    # 'boxcounting': 'counting',
    'MSD_short': 'tab:blue',
    'MSD_centre_of_mass_onepoint': 'tab:blue',
    'MSD_centre_of_mass_proximity': 'tab:blue',
    'timescaleint': 'tab:orange',
    'timescaleint_nofit': 'tab:orange',
    'boxcounting_collective': 'tab:orange',
    'boxcounting_shorttime': 'tab:orange',
    'boxcounting_first_quad': 'tab:orange',
    'C_N_simplefit': 'tab:red',
    'D0Sk': 'tab:red',
}

markers = {
    'DDM_short': '*',
    'f_short': 'x',
    'DDM_long': '*',
    'f_long': 'x',
    'DDM': '*',
    'f': '*',
    'MSD': '_',
    'MSD_short': '_',
    'MSD_long': '_',
    'MSD_first': '_',
    'MSD_centre_of_mass_onepoint': '_',
    'MSD_centre_of_mass_proximity': '_',
    'boxcounting': '_',
    'boxcounting_first_quad': '_',
    'boxcounting_shorttime': '_',
    'boxcounting_collective': 'x',
    'timescaleint': '+',
    'timescaleint_nofit': 'o',
    'D0Sk': 'o',
    'C_N_simplefit': '|',
}

def go(file, sources, PLOT_AGAINST_K, TWO_PI, logarithmic_y, fix_axes, export_destination=None):

    all_Ds = []
    
    used_sources = []

    pixel_size = None
    window_size = None

    if PRESENT_SMALL:
        figsize = (3.2, 2.8)
        if PLOT_AGAINST_K and SHOW_TWIN_K_AXIS:
            figsize = (3.2, 3.2)
    else:
        figsize = (5, 4)
        
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    MSD_D = None
    pack_frac_calced = None
    pack_frac_given  = None

    # for source in ['f', 'Fs', 'DDM', 'boxcounting', 'boxcounting_shorttime', 'MSD']:
    # for source in ['boxcounting', 'MSD', 'Fs', 'f', 'DDM']:
    # for source in ['boxcounting', f'Fs_{timescale}', f'f_{timescale}', f'DDM_{timescale}', 'MSD']:
    # for source in ['boxcounting', f'MSD_{timescale}', f'Fs_{timescale}', f'f_{timescale}', f'DDM_{timescale}', f'f_long', f'DDM_long', 'f', 'DDM', 'timescaleint']:
    # for source in [f'f_{timescale}', f'DDM_{timescale}', f'f_long', f'DDM_long', 'f', 'DDM', 'timescaleint']:
    # for source in [f'f_{timescale}', f'DDM_{timescale}', f'DDM_long', 'DDM', 'timescaleint']:
    for source in sources:

        source_label = f'{source_names[source]}'# if not source.startswith('MSD') else None
        if source in ['MSD', 'MSD_short']:
            source_label = None

        if source == 'D0Sk':
            data = common.load(f"scattering_functions/data/F_{file}.npz")
            t                 = data["t"]
            F                 = data["F"] # (num timesteps) x (num k bins)
            F_unc             = data['F_unc']
            k                 = data["k"]
            # particle_diameter = data.get('particle_diameter', 1)

            S = F[0, :]
            S_unc = F_unc[0, :]
            
            start_index = 40 # crops off k=0 delta fn

            if PLOT_AGAINST_K:
                xs = k[0, :]
            else:
                if TWO_PI:
                    xs = 2 * np.pi / k[0, :]
                    # source_label += ' $2\pi/k$'
                else:
                    xs = 1 / k[0, :]
                    # source_label += ' $1/k$'

            Ds = MSD_D / S
            D_uncs = MSD_D / S**2 * S_unc
        
        else:

            try:
                data = common.load(f'visualisation/data/Ds_from_{source}_{file}.npz')
            except FileNotFoundError:
                if source == 'MSD_short':
                    raise
                print('FileNotFound', source)
                continue
            Ds     = data['Ds']
            D_uncs = data['D_uncs']

            # assert np.all(Ds >= 0)
            # assert np.all(D_uncs >= 0)

            if Ds.size == 0:
                print(f'skipping {source} {file}, no Ds found')
                continue
            
            if source.startswith('f') or source.startswith('Fs') or source.startswith('DDM'):
                # xs = 1 / data['ks']
                # print('max k', data['ks'].max(), '2pi over over D', 2*np.pi/data['ks'].max()/2.82)

                if PLOT_AGAINST_K:
                    xs = data['ks']
                else:
                    if TWO_PI:
                        xs = 2 * np.pi / data['ks']
                        source_label += ' $2\pi/k$'
                    else:
                        xs = 1 / data['ks']
                        source_label += ' $1/k$'

                # if file == 'eleanorlong':
                #     print('REMOVING POINTS!')
                #     start = 9
                #     xs = xs[start:]
                #     Ds = Ds[start:]
                #     D_uncs = D_uncs[start:]


            elif source.startswith('boxcounting'):
                xs = data['Ls']
                if PLOT_AGAINST_K:
                    raise

                # print('THIS IS SO NAUGHTY REMOVE ME')
                # thresh = 0.576
                # Ds = Ds[xs > thresh]
                # D_uncs = D_uncs[xs > thresh]
                # xs = xs[xs > thresh]

            elif source.startswith('timescaleint') or source.startswith('C_N') :
                if PLOT_AGAINST_K:
                    print('skipping timescaleint')
                    continue
                else:
                    xs = data['Ls']

            elif source.startswith('MSD_centre_of_mass'):
                if PLOT_AGAINST_K:
                    print('skipping MSD centre of mass')
                    continue
                else:
                    xs = np.sqrt(data['Ns'] / data['density'])

            elif source.startswith('MSD'):
                MSD_D = Ds[0]
                
                # print('MSD errors hacked')
                # ax.fill_between(ax.get_xlim(), Ds[0]*0.97, Ds[0]*1.03, facecolor=colors['MSD_short'], alpha=0.5)
                # print('msd unc', D_uncs)
                # all_Ds.append(Ds[0])
                xs = np.array([])
                Ds = np.array([])
                D_uncs = np.array([])
            else:
                raise Exception(f'you need to specify the x scale for {source}')


        diameter = data.get('particle_diameter')

        
        if 'pixel_size' in data and data['pixel_size'] != None:
            if diameter and not np.isnan(diameter):
                if PLOT_AGAINST_K:
                    pixel_size = data['pixel_size'] * diameter
                else:
                    pixel_size = data['pixel_size'] / diameter
            else:
                pixel_size = data['pixel_size']
        if 'window_size_x' in data and data['window_size_x'] != None:
            if diameter and not np.isnan(diameter):
                if PLOT_AGAINST_K:
                    window_size = min(data['window_size_x'], data['window_size_y']) * diameter
                else:
                    window_size = min(data['window_size_x'], data['window_size_y']) / diameter
            else:
                window_size = min(data['window_size_x'], data['window_size_y'])

        if 'pack_frac_given' in data and data['pack_frac_given'] != None:
            pack_frac_given = data['pack_frac_given']
        if 'pack_frac_calced' in data and data['pack_frac_calced'] != None:
            pack_frac_calced = data['pack_frac_calced']

        if diameter and not np.isnan(diameter):
            if PLOT_AGAINST_K:
                xs *= diameter
            else:
                xs /= diameter
        else:
            print('not rescaling by diameter')

        assert not np.any(np.isnan(xs)), f'nan found in xs from {source}'

        # if not source.startswith('MSD'):
        color = colors[source]
        if MSD_D:
            ys = Ds/MSD_D
            yerrs = D_uncs/MSD_D
            ylabel = '$D/D_{{MSD}}$'
        else:
            ys = Ds
            yerrs = D_uncs
            ylabel = '$D$ ($\mathrm{\mu m^2/s}$)'
        lines = ax.plot(xs, ys, linestyle='none', marker=markers[source], markersize=4, color=color, label=source_label)
        ax.errorbar(xs, ys, yerr=yerrs, linestyle='none', marker='none', alpha=0.6, color=color)

        # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'
        [all_Ds.append(D) for D in ys]

        used_sources.append(source)

    ax.hlines(1, *ax.get_xlim(), color=colors['MSD_short'], linestyle='dotted', label=source_names['MSD_short'])
    print('MSD errors hacked')
    ax.fill_between(ax.get_xlim(), 0.97, 1.03, facecolor=colors['MSD_short'], alpha=0.5)

    if pack_frac_calced:
        pack_frac = pack_frac_calced
    elif pack_frac_given:
        pack_frac = pack_frac_given
    else:
        pack_frac = None
    if pack_frac:
        x = (1+pack_frac)/(1-pack_frac)**3
        ax.hlines(x, *ax.get_xlim(), label=r'$(1+\phi)/(1-\phi)^3$', color='gray')

            
    if not PLOT_AGAINST_K:
        ax.vlines(pixel_size,  min(ys), max(ys), color='gray', linestyle='dotted', label='pixel size')
        ax.vlines(window_size, min(ys), max(ys), color='gray', linestyle='dashed', label='window size')

    assert len(all_Ds) > 0, 'no Ds were found at all'

    phi = 0.34
    D_th_self = ( 1 - 1.73*phi )
    D_th_coll = ( 1 + 1.45*phi )
    # ax.hlines([D_th_coll, D_th_self], xs.min(), xs.max(), color='black')

    if logarithmic_y:
        ax.semilogy()

    ylim_expand = 1.2
    if np.nanmax(all_Ds) - np.nanmax(all_Ds) < 0.4:
        ylim_expand = 1.5
    ymin = np.nanmin(all_Ds)/ylim_expand
    if MSD_D:
        ymin = 0.3
    ax.set_ylim(ymin, np.nanmax(all_Ds)*ylim_expand)
    if fix_axes:
        ax.set_ylim(0.1, 5)
    ax.set_ylabel(ylabel)
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

    # ax.legend(fontsize=7, loc='lower right')
    ax.legend(fontsize=7)
    # ax.set_title(f'{file}, errorbars not yet all correct')
    name = '_'.join(used_sources)
    filename = f'{file}_{name}'
    if not TWO_PI:
        filename += '_oneoverk'
    
    if export_destination:
        common.save_fig(fig, f'{export_destination}/Ds_overlapped_{filename}.pdf', hide_metadata=True)
    common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_{filename}.png', dpi=200)

    print()

for file in sys.argv[1:]:

    go(file, ['MSD_short', 'D0Sk', 'f', 'f_short', 'boxcounting_shorttime', 'boxcounting_collective', 'timescaleint_nofit', 'timescaleint', 'MSD_centre_of_mass_proximity'],
        PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, fix_axes=False)

    # go(file, ['MSD_short', 'boxcounting_collective', 'timescaleint_nofit', 'timescaleint', 'C_N_simplefit'],
    #     PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, fix_axes=False)

    # go(file, ['MSD_short', 'D0Sk', 'f', 'f_short', 'f_long'],
    #     PLOT_AGAINST_K=True, TWO_PI=True, logarithmic_y=True, fix_axes=False)