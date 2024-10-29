import common
import numpy as np
import matplotlib.pyplot as plt
import sys
import countoscope_theory.structure_factor
import countoscope_theory.timescaleint
import scipy.special

from visualisation.Ds_overlapped import get_color, get_D0, get_linestyle, get_marker, source_names

SHOW_TWIN_K_AXIS = False
ERRORBAR_ALPHA = 0.5
LABELS_ON_PLOT = False
LABELS_ON_PLOT_FONTSIZE = 7
D0_SOURCE = 'MSD_first'
# D0_SOURCE = 'MSD_short'

DEFAULT_FIGSIZE = (6, 5)

def get_tau(source, file, PLOT_AGAINST_K, TWO_PI):
    data = None

    if source == 'D0Sk':
        data = common.load(f"isf/data/F_{file}.npz")
        t                 = data["t"]
        F                 = data["F"] # (num timesteps) x (num k bins)
        F_unc             = data['F_unc']
        k                 = data["k"]
        # particle_diameter = data.get('particle_diameter', 1)
        

        D_MSD, sigma, phi = get_D0(file)

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
        data = common.load(f"isf/data/F_{file}.npz")
        t                 = data["t"]
        F                 = data["F"] # (num timesteps) x (num k bins)
        F_unc             = data['F_unc']
        k                 = data["k"]
        # particle_diameter = data.get('particle_diameter', 1)
        
        D_MSD, sigma, phi = get_D0(file)
        
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

    # elif source == 'D_of_L_theory':
    #     D_MSD, sigma, phi = get_D0(file)

    #     L = np.logspace(np.log10(sigma*1e-1), np.log10(sigma*5e1), 200)
    #     L = L[::-1] # so that it's the same order as the others
    #     D = countoscope_theory.timescaleint.D_of_L(L, 1, phi, sigma)
    #     print('itssss', D[-1])

    #     xs = L
    #     Ds = D * D_MSD
    #     D_uncs = np.zeros_like(Ds)

    # elif source == 'D_of_L_theory_Lh':
    #     D_MSD, sigma, phi = get_D0(file)

    #     L = np.logspace(np.log10(sigma*1e-1), np.log10(sigma*5e1), 100)
    #     L = L[::-1] # so that it's the same order as the others
    #     D = countoscope_theory.timescaleint.D_of_L_Lh(L, 1, phi, sigma)
    #     print('itssss', D[-1])

    #     xs = L
    #     Ds = D * D_MSD
    #     D_uncs = np.zeros_like(Ds)

    elif source == 'dominiguez_theory':
        # assert '034' in file
        D_MSD, sigma, phi = get_D0(file)
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
        D_MSD, sigma, phi = get_D0(file)

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
        if source.startswith('MSD') and file == 'eleanorlong066':
            usedfile = 'eleanorlong066_div8'
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
        
        if source.startswith('f') or source.startswith('Fs') or source.startswith('DDM'):
            
            if PLOT_AGAINST_K:
                xs = data['ks']
            else:
                if TWO_PI:
                    xs = 2 * np.pi / data['ks']
                    # source_label += ' $2\pi/k$'
                else:
                    xs = 1 / data['ks']
                    # source_label += ' $1/k$'

            taus = 1 / (data['Ds'] * data['ks']**2)

        elif source.startswith('boxcounting') or source.startswith('timescaleint') or source.startswith('C_N'):
            assert False
            xs = data['Ls']
            if PLOT_AGAINST_K:
                raise

        elif source.startswith('MSD_centre_of_mass'):
            xs = np.sqrt(data['Ns'] / data['density'])

        elif source.startswith('MSD'):
            xs = []

        else:
            raise Exception(f'you need to specify the x scale for {source}')



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
                    window_size = min(data['window_size_x'], data['window_size_y']) * diameter
                else:
                    window_size = min(data['window_size_x'], data['window_size_y']) / diameter
            else:
                window_size = min(data['window_size_x'], data['window_size_y'])

        if 'pack_frac_given' in data and data['pack_frac_given'] != None:
            pack_frac_given = data['pack_frac_given']
        if 'pack_frac_calced' in data and data['pack_frac_calced'] != None:
            pack_frac_calced = data['pack_frac_calced']

    assert not np.any(np.isnan(xs)), f'nan found in xs from {source}'
    
    return xs, taus, np.zeros_like(taus), pixel_size, window_size, pack_frac_given, pack_frac_calced, diameter


def go(file, sources, PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, ylim=None,
       output_filename=None, export_destination=None, show_pixel=True,
       show_window=True, show_pack_frac_plateau=False, figsize=DEFAULT_FIGSIZE,
       label_k_scaling=False, label_pack_frac=False, hide_msd=False, ax=None):

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
    taus = {}
    tau_uncs = {}

    used_sources = []

    # get all the data
    for source in sources:
        try:
            xs[source], taus[source], tau_uncs[source], pixel_size_temp, window_size_temp, pack_frac_given_temp, pack_frac_calced_temp, diameter_temp = get_tau(source, file, PLOT_AGAINST_K, TWO_PI)
        
            if pixel_size_temp:       pixel_size       = pixel_size_temp
            if window_size_temp:      window_size      = window_size_temp
            if pack_frac_calced_temp: pack_frac_calced = pack_frac_calced_temp
            if pack_frac_given_temp:  pack_frac_given  = pack_frac_given_temp
            if diameter_temp:         diameter         = diameter_temp

            used_sources.append(source)
            found_at_least_one_file = True

        except FileNotFoundError:
            print('FileNotFound', source)

    assert found_at_least_one_file

    if D0_SOURCE in used_sources:
        # ax.set_ylabel('$D/D_{{MSD}}$')
        rescale_y = taus[D0_SOURCE][0]
    else:
        ax.set_ylabel(r'$\tau$ ($\mathrm{s}$)')
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
            D_MSD, _, _ = get_D0(file)
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

        try:
            print('x', source, taus[source][np.nanargmin(xs[source])])
        except:
            pass
        print('X', source, taus[source][0])
        xs[source] /= rescale_x
        # print('rescale_y', rescale_y)

        color = get_color(source)
        ys = taus[source] / rescale_y
        yerrs = tau_uncs[source] / rescale_y
        
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
            print('ys', source, ys)
            ax.plot(xs[source], ys, linestyle=get_linestyle(source), marker=get_marker(source), markersize=4, color=color,
                label=plot_label, zorder=zorder)
            ax.errorbar(xs[source], ys, yerr=yerrs, linestyle='none', marker='none', alpha=ERRORBAR_ALPHA, color=color)

        if LABELS_ON_PLOT:
                
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



        # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'
        [all_Ds.append(D) for D in ys]


    if not PLOT_AGAINST_K:
        if show_pixel:
            ax.vlines(pixel_size,  min(ys), max(ys), color='gray', linestyle='dotted', label='pixel size')
        if show_window:
            ax.vlines(window_size, min(ys), max(ys), color='gray', linestyle='dashed', label='window size')

    assert len(all_Ds) > 0, 'no Ds were found at all'

    # get max time
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    time_step = data['time_step']
    max_time = particles[:, 2].max() * time_step
    print(max_time, 'max_time', ax.get_xlim())
    ax.hlines(max_time, *ax.get_xlim(), label='max time')

    if logarithmic_y:
        ax.semilogy()

    # ylim_expand = 1.2
    # if np.nanmax(all_Ds) - np.nanmax(all_Ds) < 0.4:
    #     ylim_expand = 1.5
    # ymin = np.nanquantile(all_Ds, 0.02)/ylim_expand
    # ax.set_ylim(ymin, np.nanquantile(all_Ds, 0.98)*ylim_expand)
    
    ax.semilogx()
    if PLOT_AGAINST_K:
        

        if SHOW_TWIN_K_AXIS:
            realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
            # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
            realspace_ax.set_xlabel(r'$2\pi/k / \sigma$')
    else:
        ax.set_xlabel(r'$L / \sigma$')
        
    if not LABELS_ON_PLOT:
        ax.legend(fontsize=5)

    if label_k_scaling:
        k_scaling = r'$k \rightarrow L = 2\pi/k$' if TWO_PI else r'$k \rightarrow L = 1/k$'
        ax.text(0.12, 0.8, k_scaling, transform=ax.transAxes)

    # if file == 'eleanorlong066':
    #     D0 = 0.015
    #     ax.hlines(D0, *ax.get_xlim())
    #     s = (1+0.66)/(1-0.66)**3
    #     ax.hlines(D0*s, *ax.get_xlim())

    
    if not ax_supplied:
        if export_destination:
            common.save_fig(fig, export_destination, hide_metadata=True)
        if output_filename:
            common.save_fig(fig, f'visualisation/figures_png/taus_overlapped_{output_filename}.png', dpi=200)

    print()


if __name__ == '__main__':
    for file in sys.argv[1:]:

        filename = file

        go(file, [
            # 'MSD_first',
                #   'MSD_long',
                #   'D0Sk', 
                #   'MSD_short',
                #   'D0Sk_theory',
                #   'dominiguez_theory',
                #   'panzuela_theory',
                #   'f_D1', 'f_D2',
                  'f', 
                'f_short',
                'f_long', 
                  'f_first',
                # 'f_short',
                #   'boxcounting_shorttime',
                    # 'boxcounting_first_quad',
                    # 'boxcounting_collective',
                #   'MSD_centre_of_mass_proximity',
                # 'D_of_L_theory',
                # 'D_of_L_theory_Lh',
            ],
            PLOT_AGAINST_K=True, TWO_PI=True, logarithmic_y=True, output_filename=filename)
