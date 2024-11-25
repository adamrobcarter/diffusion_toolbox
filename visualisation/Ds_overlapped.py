import common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys
import countoscope_theory.structure_factor
import countoscope_theory.timescaleint
import scipy.special

SHOW_TWIN_K_AXIS = False
ERRORBAR_ALPHA = 0.4
LABELS_ON_PLOT = False
LABELS_ON_PLOT_FONTSIZE = 7
LEGEND_FONTSIZE = 7
D0_SOURCE = 'MSD_first'
# D0_SOURCE = 'MSD_short'
PREVENT_Y_RESCALING = False
YLIM = None
DEFAULT_FIGSIZE = (6, 5)

source_names = {
    'DDM': 'DDM',
    'DDM_short': 'DDM short',
    'DDM_long': 'DDM long',
    'f': '$f(k, \Delta t)$',
    'Fs': '$F_s(k, \Delta t)$',
    'f_short': '$f(k, \Delta t)$ short',
    'Fs_short': '$F_s(k, \mathrm{short})$',
    'f_long': '$f(k, \Delta t)$ long',
    'f_first': '$f(k, \Delta t)$ first point',
    'Fs_long': '$F_s(k, \mathrm{long})$',
    'boxcounting': 'counting full fit',
    'MSD': 'MSD',
    'MSD_short': 'MSD',
    'MSD_long': 'MSD long time',
    'MSD_first': 'MSD first',
    'MSD_centre_of_mass_onepoint': 'MSD CoM',
    'MSD_centre_of_mass_proximity': 'MSD CoM prox',
    'boxcounting_shorttime': 'Countoscope short time fit',
    'boxcounting_first_quad': 'Countoscope first point',
    'boxcounting_collective': 'Countoscope full fit',
    'timescaleint': 'timescale integral',
    'timescaleint_nofit': 'timescale integral (no fit)',
    'D0Sk': r'$D_{MSD}/S(k)$',
    'D0Sk_theory': r'$D(k)$ theory (no hydro)',
    'C_N_simplefit': '$C_N$ fit',
    'D_of_L_theory': '$D(L)$ theory (no hydro)',
}

colors = {
    'DDM_short': 'tab:purple',
    'f': 'lime',
    'f_short': 'tab:green',
    'f_long': 'yellowgreen',
    'f_first': 'green',
    'Fs_short': 'tab:green',
    'f_D1': 'tab:green',
    'f_D2': 'turquoise',
    # 'boxcounting': 'counting',
    'MSD_first': 'tab:cyan',
    'MSD_short': 'tab:blue',
    'MSD_centre_of_mass_onepoint': 'tab:blue',
    'MSD_centre_of_mass_proximity': 'tab:blue',
    'timescaleint': 'gold',
    'timescaleint_nofit': 'tab:orange',
    'timescaleint_nofit_cropped': 'tab:orange',
    'boxcounting_collective': 'tab:orange',
    'boxcounting_shorttime': 'tab:orange',
    'boxcounting_first_quad': 'tab:orange',
    'C_N_simplefit': 'tab:red',
    'D0Sk': 'mediumspringgreen',
    'D0Sk_theory': 'aquamarine', # tab:olive,
    'D_of_L_theory': 'tab:red',
    'NtN0_first': 'deeppink',
}

markers = {
    'DDM_short': '*',
    'f_short': 'x',
    'DDM_long': '*',
    'f_long': 'x',
    'DDM': '*',
    'f': '*',
    'f_first': '*',
    'MSD': '_',
    'MSD_short': '_',
    'MSD_long': '_',
    'MSD_first': '_',
    'MSD_centre_of_mass_onepoint': '_',
    'MSD_centre_of_mass_proximity': '_',
    'boxcounting': '_',
    'boxcounting_first_quad': '_',
    'boxcounting_shorttime': '*',
    'boxcounting_collective': 'x',
    'timescaleint': 'o',
    'timescaleint_nofit': 'o',
    'timescaleint_nofit_cropped': 'o',
    'D0Sk': 'o',
    'D0Sk_theory': 'none',
    'C_N_simplefit': '|',
    'D_of_L_theory': 'none',
    'dominiguez_theory': 'none',
    'D_of_L_theory_Lh': 'none',
}

linestyle = {
    'D0Sk_theory': '-',
    'D_of_L_theory': '-',
    'D_of_L_theory_Lh': '-',
    'dominiguez_theory': '-',
}


def get_D0(file):
    # we don't do these ones in get_D0_filename cause sometimes we might want to compare these ones
    suffixes = ['_crop', '_trim', '_first', '_smallbins', '_nozero', '_no_overlap', '_long', '_longer']
    for suffix in suffixes:
        if suffix in file:
            file = file.split(suffix)[0]
            
    file = get_D0_filename(file)
            
    data = common.load(f"visualisation/data/Ds_from_{D0_SOURCE}_{file}.npz")
    D_MSD = data["Ds"][0]
    sigma = data['particle_diameter']
    if 'pack_frac' in data:
        phi = data['pack_frac']
    else:
        print('using pack frac *given*')
        phi = data['pack_frac_given']

    print(f'D_MSD = {D_MSD}')
    return D_MSD, sigma, phi

def get_D0_filename(file):

    suffixes = ['_merged'] # these are ones we will never want to compare
    for suffix in suffixes:
        if suffix in file:
            file = file.split(suffix)[0]

    if file in ['eleanorlong066', 'sim_nohydro_034_L640', 'sim_nohydro_034_L1280',
                'sim_nohydro_010_L320', 'sim_nohydro_010_L544', 'sim_nohydro_010_L544_dt2', 'sim_nohydro_010_L640',
                'brennan_hydro_010_L544', 'brennan_hydro_010_L1280',
                'sim_nohydro_011_L320', 'sim_nohydro_011_L320_long',]:
        file += '_div8'
    return file

def get_L_and_D(source, file, PLOT_AGAINST_K, TWO_PI, D_MSD, phi, sigma):
    data = None

    if source == 'D0Sk':
        data = common.load(f"isf/data/F_{file}.npz")
        t                 = data["t"]
        F                 = data["F"] # (num timesteps) x (num k bins)
        F_unc             = data['F_unc']
        k                 = data["k"]

        S = F[0, :]
        S_unc = F_unc[0, :]
        
        # min_index = np.nanargmin(S)

        if PLOT_AGAINST_K:
            xs = k[0, :]
        else:
            if TWO_PI:
                xs = 2 * np.pi / k[0, :]
                # source_label += ' $2\pi/k$'
            else:
                xs = 1 / k[0, :]
                # source_label += ' $1/k$'

        Ds = D_MSD / S[:]
        D_uncs = D_MSD / S[:]**2 * S_unc[:]

    elif source == 'D0Sk_theory':
        
        L = np.logspace(np.log10(sigma*1e-1), np.log10(sigma*3e1), 100)
        L = L[::-1] # so that it's the same order as the others
        k = 2*np.pi/L
        S = countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)

        if PLOT_AGAINST_K:
            xs = k
        else:
            if TWO_PI:
                xs = 2 * np.pi / k
                # source_label += ' $2\pi/k$'
            else:
                xs = 1 / k
                # source_label += ' $1/k$'

        Ds = D_MSD / S
        D_uncs = np.zeros_like(Ds)

    elif source == 'D_of_L_theory':

        L = np.logspace(np.log10(sigma*0.5e-1), np.log10(sigma*5e1), 100)
        L = L[::-1] # so that it's the same order as the others
        D = countoscope_theory.timescaleint.D_of_L(L, D_MSD, phi, sigma)

        xs = L
        Ds = D
        D_uncs = np.zeros_like(Ds)

    elif source == 'D_of_L_theory_Lh':

        L = np.logspace(np.log10(sigma*0.5e-1), np.log10(sigma*5e1), 100)
        L = L[::-1] # so that it's the same order as the others
        D = countoscope_theory.timescaleint.D_of_L_Lh(L, 1, phi, sigma)

        xs = L
        Ds = D * D_MSD
        D_uncs = np.zeros_like(Ds)

    elif source == 'dominiguez_theory':
        # xs = np.logspace()
        k = np.logspace(-2, 0.5)
        # L = 2*np.pi/k
        L = 1/k
        print('used k=1/L')
        # L*= 1.5
        xs = L
        Ds = np.zeros_like(xs)
        Lh = sigma / (3 * phi) # Dominiguez 2014 eq 31a
        # Lh = 0.2 * sigma
        # Lh = sigma
        Ds[L <= Lh] = D_MSD
        Ds[L > Lh] = D_MSD / (k[L > Lh] * Lh)
        D_uncs = np.zeros_like(Ds)

    elif source == 'panzuela_theory':

        L = np.logspace(-1, 2, 200)
        k = 1/L
        delta = 0.08 # Gaussian trap width
        print('used k=1/L')
        H = 3 * phi * delta / sigma * (1 / (2*k*delta + k*delta) * np.exp(k**2 * delta**2) * scipy.special.erfc(k*delta) - 1/np.sqrt(np.pi))
        S = countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)
        Ds = D_MSD * H / S
        D_uncs = np.zeros_like(Ds)
        xs = L

    else:
        usedfile = file
        if source.startswith('MSD'):
            usedfile = get_D0_filename(file)
        data = common.load(f'visualisation/data/Ds_from_{source}_{usedfile}.npz')


        Ds     = data['Ds']
        D_uncs = data['D_uncs']

        if source.startswith('MSD'):
            print('loaded', source, Ds)

        # assert np.all(Ds >= 0)
        # assert np.all(D_uncs >= 0)

        if Ds.size == 0:
            print(f'skipping {source} {file}, no Ds found')
            raise FileNotFoundError(f'skipping {source} {file}, no Ds found')
        
        if source.startswith('f') or source.startswith('F_s') or source.startswith('DDM'):
            
            if PLOT_AGAINST_K:
                xs = data['ks']
            else:
                if TWO_PI:
                    xs = 2 * np.pi / data['ks']
                    # source_label += ' $2\pi/k$'
                else:
                    xs = 1 / data['ks']
                    # source_label += ' $1/k$'
            print(source, 'min f', data['ks'].min(), '2pi/min f', 2*np.pi/data['ks'].min())

        elif source.startswith('boxcounting') or source.startswith('timescaleint') or source.startswith('C_N') or source.startswith('NtN0'):
            xs = data['Ls']
            T_of_L = countoscope_theory.timescaleint.T_of_L(xs, D_MSD, phi, sigma)
            print('L', xs)
            print('T', T_of_L)
            max_time = data['max_time_hours'] * 60 * 60
            print('max time', max_time)
            print(10 * T_of_L < max_time)

            num_before = xs.size
            xs     = xs    [   10 * T_of_L < max_time]
            Ds     = Ds    [   10 * T_of_L < max_time]
            print(D_uncs.shape)
            if len(D_uncs.shape) == 1:
                D_uncs = D_uncs[10 * T_of_L < max_time]
            else:
                D_uncs = D_uncs[:, 10 * T_of_L < max_time]
            assert xs.size, 'we dropped all points'
            print(xs)
            print(f'kept L<{xs[-1]/sigma:.2g}σ ({xs.size/num_before:.2f})')
            
            assert not PLOT_AGAINST_K

        elif source.startswith('MSD_centre_of_mass'):
            xs = np.sqrt(data['Ns'] / data['density'])

        elif source.startswith('MSD'):
            xs = []

        else:
            raise Exception(f'you need to specify the x scale for {source}')

    print(f'avg D_unc/D0 = {np.nanmean(D_uncs/Ds):.3f}, max =  {np.nanmax(D_uncs/Ds):.3f}')

    pixel_size = None
    window_size = None
    pack_frac_calced = None
    pack_frac_given = None
    diameter = None
    
    if data:
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
                    window_size = max(data['window_size_x'], data['window_size_y']) * diameter
                else:
                    window_size = max(data['window_size_x'], data['window_size_y']) / diameter
            else:
                window_size = max(data['window_size_x'], data['window_size_y'])

        if 'pack_frac_given' in data and data['pack_frac_given'] != None:
            pack_frac_given = data['pack_frac_given']
        if 'pack_frac_calced' in data and data['pack_frac_calced'] != None:
            pack_frac_calced = data['pack_frac_calced']

    assert not np.any(np.isnan(xs)), f'nan found in xs from {source}'
    
    if np.all(Ds == 0):
        print('all Ds zero!')
    if np.all(D_uncs == 0):
        print('all D_uncs zero!')

    return xs, Ds, D_uncs, pixel_size, window_size, pack_frac_given, pack_frac_calced, diameter


def go(file, sources, PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, ylim=None,
       output_filename=None, export_destination=None, show_pixel=True,
       show_window=True, show_pack_frac_plateau=False, figsize=DEFAULT_FIGSIZE,
       label_k_scaling=False, label_pack_frac=False, hide_msd=False, ax=None,
       legend_fontsize=None):

    all_Ds = []
    found_at_least_one_file = False
        
    ax_supplied = ax != None
    if not ax_supplied:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    pack_frac_calced = None
    pack_frac_given  = None
    pixel_size = None
    window_size = None
    diameter = None

    xs = {}
    Ds = {}
    D_uncs = {}

    used_sources = []

    D_MSD, sigma, phi = get_D0(file)

    # get all the data
    for source in sources:
        try:
            xs[source], Ds[source], D_uncs[source], pixel_size_temp, window_size_temp, \
                pack_frac_given_temp, pack_frac_calced_temp, diameter_temp = \
                get_L_and_D(source, file, PLOT_AGAINST_K, TWO_PI, D_MSD, phi, sigma)
        
            if pixel_size_temp:       pixel_size       = pixel_size_temp
            if window_size_temp:      window_size      = window_size_temp
            if pack_frac_calced_temp: pack_frac_calced = pack_frac_calced_temp
            if pack_frac_given_temp:  pack_frac_given  = pack_frac_given_temp
            if diameter_temp:         diameter         = diameter_temp

            used_sources.append(source)
            found_at_least_one_file = True

        except FileNotFoundError as e:
            print('FileNotFound', e)

    assert found_at_least_one_file

    if D0_SOURCE in used_sources and not PREVENT_Y_RESCALING:
        ax.set_ylabel(r'$D/D_\mathrm{self}$')
        rescale_y = Ds[D0_SOURCE][0]
    else:
        ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
        rescale_y = 1

    if diameter and not np.isnan(diameter):
        if PLOT_AGAINST_K:
            rescale_x = 1/diameter
            ax.set_xlabel(r'$k \sigma$')
        else:
            rescale_x = diameter
            ax.set_xlabel(r'$L/\sigma$')
    else:
        rescale_x = 1
        if PLOT_AGAINST_K:
            ax.set_xlabel(r'$k$')
        else:
            ax.set_xlabel(r'$L$')

    assert sum([len(x) for x in xs.values()]) > 0, 'xs was empty'
    xmin = min([min(x) for x in xs.values() if len(x)>0]) / rescale_x
    xmax = max([max(x) for x in xs.values() if len(x)>0]) / rescale_x

    if pack_frac_calced:
        pack_frac = pack_frac_calced
    elif pack_frac_given:
        pack_frac = pack_frac_given
    else:
        pack_frac = None
    if pack_frac:
        if show_pack_frac_plateau:
            x = (1+pack_frac)/(1-pack_frac)**3 * D_MSD
            Dc_theory_color = 'gray'
            Dc_theory_label = r'$(1+\phi)/(1-\phi)^3$' if not LABELS_ON_PLOT else None
            ax.hlines(x, xmin, xmax, label=Dc_theory_label, color=Dc_theory_color, linestyle=':')
            if LABELS_ON_PLOT:
                ax.text(xmax, x/1.1, '$D_c$ theory (no hydro)', ha='right', va='top', color=Dc_theory_color, fontsize=LABELS_ON_PLOT_FONTSIZE)

        if label_pack_frac:
            ax.text(0.12, 0.8, f'$\phi={pack_frac:.2f}$', transform=ax.transAxes, color=common.FIT_COLOR, fontsize=10)

    # do the actual plotting
    for source in used_sources:
        source_label = f'{source_names.get(source, source)}'# if not source.startswith('MSD') else None

        xs[source] = xs[source] / rescale_x

        color = get_color(source)
        ys = Ds[source] / rescale_y
        yerrs = D_uncs[source] / rescale_y
        
        assert not np.any(np.isnan(xs[source]))

        plot_label = source_label if not LABELS_ON_PLOT else None
        if source == 'D_of_L_theory':
            print('PLOTTING')
            # print(xs[source], ys, linestyle.get(source, 'none'))
        zorder = -1 if source.endswith('theory') else None

        if source.startswith('MSD'):
            if not hide_msd:
                ax.hlines(ys[0], xmin, xmax, color=color, linestyle='dotted', label=plot_label)
                print('MSD errors hacked')
                ax.fill_between(ax.get_xlim(), ys[0]*0.97, ys[0]*1.03, facecolor=color, alpha=0.5)

        else:
            # print('ys', source, ys)
            print('for', file, source, 'plotting', xs[source].shape)
            print('getting marker', file, source, get_marker(source))
            ax.plot(xs[source], ys, linestyle=get_linestyle(source), marker=get_marker(source), markersize=4, color=color,
                label=plot_label, zorder=zorder)
            ax.errorbar(xs[source], ys, yerr=yerrs, linestyle='none', marker='none', alpha=ERRORBAR_ALPHA, color=color)

        if LABELS_ON_PLOT:
                
            # old code for MSD labels:
            #         xpos = xmin if PLOT_AGAINST_K else xmax
            #         ha = 'left' if PLOT_AGAINST_K else 'right'
            #         ax.text(xpos, ys[0]/1.2, source_label, ha=ha, va='top', color=colors[source], fontsize=LABELS_ON_PLOT_FONTSIZE)

            offset = len(ys) // 6 if PLOT_AGAINST_K else -len(ys) // 6
            if source == 'D_of_L_theory':
                offset = int(len(ys) * 0.6)
            if source == 'D0Sk_theory':
                offset = len(ys) // 2
            if source == 'boxcounting_collective':
                offset = -len(ys) // 3
            xpos = xs[source][offset] * 1.2
            ypos = ys[offset]*1.2
            ha = 'left' if PLOT_AGAINST_K else 'right'
            if source == 'f_short' and PLOT_AGAINST_K == False:
                xpos = xs[source][5] * 1.5
                ypos = ys[5]*1.2
                ha = 'left'
            ax.text(xpos, ypos, source_label, ha=ha, va='bottom', color=color, fontsize=LABELS_ON_PLOT_FONTSIZE)

        if np.all(yerrs == 0):
            print(source, 'all errors zero')


        print('Ds', file, source, ys)
        # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'
        [all_Ds.append(D) for D in ys]



    assert len(all_Ds) > 0, 'no Ds were found at all'

    # phi = 0.34
    # D_th_self = ( 1 - 1.73*phi )
    # D_th_coll = ( 1 + 1.45*phi )

    if logarithmic_y:
        ax.semilogy()
        ax.yaxis.set_minor_formatter(matplotlib.ticker.LogFormatter()) # prevent scientific notation on axes
        ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatter()) # prevent scientific notation on axes

    if not ylim:
        ylim_expand = 1.2
        if np.nanmax(all_Ds) - np.nanmax(all_Ds) < 0.4:
            ylim_expand = 1.5
        ylim = (np.nanquantile(all_Ds, 0.05)/ylim_expand, np.nanquantile(all_Ds, 0.95)*ylim_expand)
        # if ylim:
        #     ax.set_ylim(*ylim)
        # ax.set_ylabel('$D/D_0$')
        # ax.set_xticks([])
    ax.set_ylim(*ylim)

    if not PLOT_AGAINST_K:
        if show_pixel:
            ax.vlines(pixel_size,  *ylim, color='gray', linestyle='dotted', label='pixel size')
        if show_window:
            ax.vlines(window_size, *ylim, color='gray', linestyle='dashed', label='window size (max)')



    ax.semilogx()
    if PLOT_AGAINST_K:
        

        if SHOW_TWIN_K_AXIS:
            realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
            # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
            realspace_ax.set_xlabel(r'$2\pi/k / \sigma$')
    else:
        ax.set_xlabel(r'$L / \sigma$')
    # ax.set_xlabel(r'$k / (2\pi / \sigma)$')
    # ax.semilogy()

    if not LABELS_ON_PLOT:
        ax.legend(fontsize=legend_fontsize)
    # ax.set_title(f'{file}, errorbars not yet all correct')
    # name = '_'.join(used_sources)
    # filename = f'{file}_{name}'
    # if not TWO_PI:
    #     filename += '_oneoverk'

    if label_k_scaling:
        k_scaling = r'$k \rightarrow L = 2\pi/k$' if TWO_PI else r'$k \rightarrow L = 1/k$'
        ax.text(0.12, 0.8, k_scaling, transform=ax.transAxes)
    
    if not ax_supplied:
        if export_destination:
            common.save_fig(fig, export_destination, hide_metadata=True)
        # if output_filename:
        common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_{output_filename}.png', dpi=200)

    print()

def get_linestyle(source):
    return linestyle.get(trim_source(source), 'none')

def get_marker(source):
    print('marker', source, trim_source(source), markers.get(trim_source(source)))
    return markers.get(trim_source(source), 'o')

def get_color(source):
    return colors.get(trim_source(source), common.FIT_COLOR)

suffixes = ['_obs', '_nmsdfitinter', '_nmsdfit', '_var', '_varmod', '_sDFT']

def trim_source(source):
    for suffix in suffixes:
        source = source.split(suffix)[0]
    return source

if __name__ == '__main__':
    for file in sys.argv[1:]:

        filename = file

        go(file, ['MSD_first',
                #   'MSD_long',
                #   'D0Sk', 
                #   'MSD_short',
                  'D0Sk_theory',
                #   'dominiguez_theory',
                #   'panzuela_theory',
                #   'f_D1', 'f_D2',
                #   'f', 'f_short', 'f_long', 
                #   'f_first',
                # 'f_short',
                'f_long',
                'f',
                'f_first_first',
                # 'F_s',
                'F_s_first',
                # 'F_s_long',
                #   'boxcounting_shorttime',
                #     'boxcounting_first_quad',
                    'boxcounting_collective',
                #   'timescaleint_nofit_cropped_var',
                #   'timescaleint_var',
                  'timescaleint_nmsdfitinter',
                #   'MSD_centre_of_mass_proximity',
                'D_of_L_theory',
                # 'D_of_L_theory_Lh',
                # 'NtN0_first',
                # 'NtN0_fit',
                # 'F_s_first',
            ],
            PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, output_filename=filename,
            ylim=YLIM, legend_fontsize=LEGEND_FONTSIZE
            )

        # go(file, ['MSD_short', 'boxcounting_collective', 'timescaleint_nofit', 'timescaleint', 'C_N_simplefit'],
        #     PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True)

        # go(file, ['MSD_short', 'D0Sk', 'f', 'f_short', 'f_long'],
        #     PLOT_AGAINST_K=True, TWO_PI=True, logarithmic_y=True)