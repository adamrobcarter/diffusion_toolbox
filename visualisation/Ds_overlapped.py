import common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys, warnings, os
import countoscope_theory.structure_factor
import countoscope_theory.timescaleint
import scipy.special

import visualisation.Ds_overlapped_mult

D0_SOURCE = 'MSD_first'
# D0_SOURCE = 'MSD_short'

if __name__ == '__main__':
    for file in sys.argv[1:]:
        fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.5))

        sources = [
                # 'MSD_first',
                #   'MSD_long',
                #   'D0Sk', 
                #   'MSD_short',
                  'D0Sk_theory',
                #   'DDM',
                #   'dominiguez_theory',
                #   'panzuela_theory',
                #   'f_D1', 'f_D2',
                #   'f', 'f_short', 'f_long', 
                #   'f_first',
                # 'f_short',
                # 'f_long',
                # 'f',
                # 'f_first',
                # 'f_first_first',
                # 'f_t0.5',
                # 'f_t2',
                # 'f_t8',
                # 'f_t16',
                # 'f_t32',
                'f_t64',
                'f_t256',
                'f_t1024',
                'f_t4092',
                # 'f_t16384',
                # 'F_first32_first',
                # 'F_s',
                # 'F_s_first',
                # 'F_s_long',
                #   'boxcounting_shorttime',
                    # 'boxcounting_first_quad',
                    # 'boxcounting_collective_var',
                #   'timescaleint_',
                #   'timescaleint_nmsdfitinter',

                  'timescaleint_nofit_cropped_var',

                #   'timescaleint_nofit_cropped_sDFT',
                #   'timescaleint_fixexponent_cutoff',
                # 'timescaleint_fixexponent_var',
                # 'timescaleint_fixexponent_sDFT',
                # 'timescaleint_fixexponent_target_fixexponent',
                #   'timescaleint_nofit_cropped_noshort_var',
                #   'timescaleint_nofit_cropped_nmsdfitinter',
                #   'MSD_centre_of_mass_proximity',
                # 'D_of_L_theory',
                # 'D_of_L_theory_Lh',
                # 'NtN0_first',
                # 'NtN0_fit',
                # 'F_s_first',
            ]

        colors = [[common.colormap(i/len(sources)) for i in range(len(sources))]]
        print(colors)

        visualisation.Ds_overlapped_mult.go(
            colors  = colors,
            ax      = ax,
            files   = [file], 
            sources = sources,
            # plot_against_k=False,
            # logarithmic_y=True, output_filename=filename,
            # ylim=YLIM,
            legend_fontsize=8,
            linestyles = ['-']
            )
        ax.set_ylim(0.7, 6)

        # common.add_exponential_index_indicator(ax, 1/3, (8, 2), 'L')    
        
        common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_{file}.png')

        # go(file, ['MSD_short', 'boxcounting_collective', 'timescaleint_nofit', 'timescaleint', 'C_N_simplefit'],
        #     PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True)

        # go(file, ['MSD_short', 'D0Sk', 'f', 'f_short', 'f_long'],
        #     PLOT_AGAINST_K=True, TWO_PI=True, logarithmic_y=True)



"""

SHOW_TWIN_K_AXIS = False
ERRORBAR_ALPHA = 0.4
LABELS_ON_PLOT = False
LABELS_ON_PLOT_FONTSIZE = 7
LEGEND_FONTSIZE = 7
PREVENT_Y_RESCALING = False
YLIM = None
DEFAULT_FIGSIZE = (6, 5)
"""
source_names = {
    'DDM': 'DDM',
    'DDM_short': 'DDM short',
    'DDM_long': 'DDM long',
    'f': '$f(k, t)$',
    'Fs': '$F_s(k, t)$',
    'f_short': '$f(k, t)$ short',
    'Fs_short': '$F_s(k, \mathrm{short})$',
    'f_long': '$f(k, t)$ long',
    'f_first': '$f(k, t)$ first point',
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
    'f_first_first': 'lime',
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
    suffixes = ['_crop', '_trim', '_first', '_smallbins', '_nozero', '_no_overlap',
                '_long', '_longer', '_moreoverlap', '_spacing', '_frac_of_window', '_windowed', '_nowindow', '_bhwindow',
                # '_mixt'
                '_xk',# '_unmix'
                ]
    if '_mixt' in file:
        warnings.warn('allowing mixt for now, would be good to do a proper msd calculation though')
    for suffix in suffixes:
        if suffix in file:
            file = file.split(suffix)[0]

    if '_unmix' in file:
        file = file.split('_unmix')[0] + '_mixt'
            
    file = get_D0_filename(file)
            
    data = common.load(f"visualisation/data/Ds_from_{D0_SOURCE}_{file}.npz")
    D_MSD = data["Ds"][0]
    sigma = data['particle_diameter']
    if 'pack_frac' in data:
        phi = data['pack_frac']
    else:
        print('using pack frac *given*')
        if 'pack_frac_given' in data:
            phi = data['pack_frac_given']
        else:
            warnings.warn('pack frac given not in data')
            phi = np.nan

    print(f'D_MSD = {D_MSD}')
    return D_MSD, sigma, phi

def get_D0_filename(file):
    print(f'get_D0_filename({file})')

    if file in ['sim_nohydro_011_L320_test_singlet_mixt', 'sim_nohydro_011_L320_test_mixt', 'sim_nohydro_011_L320_test_singlet']: # remove after this problem solved
        return 'sim_nohydro_011_L640_div64'

    suffixes = ['_merged', '_crop1.0'] # these are ones we will never want to compare
    for suffix in suffixes:
        if suffix in file:
            file = file.split(suffix)[0]

    root = 'visualisation/data/Ds_from_MSD_first'

    if os.path.isfile(f'{root}_{file}.npz'):
        return file
    elif os.path.isfile(f'{root}_{file}_div8.npz'):
        return f'{file}_div8'
    elif os.path.isfile(f'{root}_{file}_div64.npz'):
        return f'{file}_div64'
    else:
        raise Exception(f'MSD file not found for {file} ({root}_{file}_div64.npz not found)')

    if file in ['eleanorlong066', 'sim_nohydro_034_L640', 'sim_nohydro_034_L1280',
                'sim_nohydro_010_L320', 'sim_nohydro_010_L544', 'sim_nohydro_010_L544_dt2', 'sim_nohydro_010_L640',
                'brennan_hydro_010_L544', 'brennan_hydro_010_L1280',
                'sim_nohydro_011_L320', 'sim_nohydro_011_L320_long', 'sim_nohydro_011_L320_longer',
                'sim_nohydro_002_L320', 'sim_nohydro_002_L320_long', 'sim_nohydro_002_L320_longer',
                'sim_nohydro_002_L640', 'sim_nohydro_011_L640',
                'sim_nohydro_011_L160',
                'brennan_hydro_002_L320', 'brennan_hydro_002_L640', 'brennan_hydro_011_L320']:
        file += '_div8'

    if file in ['sim_nohydro_011_L1280']:
        file += '_div64'

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
        
        L = np.logspace(np.log10(sigma*1e-1), np.log10(sigma*100), 100)
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

        if Ds.size == 0:
            print(f'skipping {source} {file}, no Ds found')
            raise FileNotFoundError(f'skipping {source} {file}, no Ds found')
        
        if source.startswith('f') or source.startswith('F_s') or source.startswith('F_first') or source.startswith('DDM'):
            
            if PLOT_AGAINST_K:
                xs = data['ks']
            else:
                if TWO_PI:
                    xs = 2 * np.pi / data['ks']
                    # source_label += ' $2\pi/k$'
                else:
                    xs = 1 / data['ks']
                    # source_label += ' $1/k$'

        elif source.startswith('boxcounting') or source.startswith('timescaleint') or source.startswith('C_N') or source.startswith('NtN0'):
            xs = data['Ls']
            # max_time = data['max_time_hours'] * 60 * 60

            # max_L = np.sqrt(4 * D_MSD * max_time / 10)
            # keep = xs < max_L
            # T_of_L = countoscope_theory.timescaleint.T_of_L(xs, D_MSD, phi, sigma)
            # keep = 10 * T_of_L < max_time

            num_before = xs.size
            # xs = xs[keep]
            # Ds = Ds[keep]
            # if len(D_uncs.shape) == 1:
            #     D_uncs = D_uncs[keep]
            # else:
            #     D_uncs = D_uncs[:, keep]
            # assert xs.size, 'we dropped all points'
            
            print(f'kept L<{xs[-1]/sigma:.2g}σ ({xs.size/num_before:.2f})')
            
            assert not PLOT_AGAINST_K

        elif source.startswith('MSD_centre_of_mass'):
            Ls = np.sqrt(data['Ns'] / data['density'])
            ks = 2*np.pi/Ls
            Ds = Ds / data['Ns']
            Ds = Ds / common.structure_factor_2d_hard_spheres(ks, 0.34, 3)
            # Ds = D_MSD / Ds
            print('Ds', Ds)

            xs = Ls

        elif source.startswith('MSD'):
            xs = []

        else:
            raise Exception(f'you need to specify the x scale for {source}')

    print(f'avg D_unc/D0 = {np.nanmean(D_uncs/Ds):.3f}, max =  {np.nanmax(D_uncs/Ds):.3f}')

    pixel_size = None
    window_size = None
    pack_frac_calced = None
    pack_frac_given = None
    
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

    if not ('theory' in source) and False:
        print('forcing min 3%% errors')
        errs_too_low = D_uncs/np.abs(Ds) < 0.03
        if len(D_uncs.shape) == 1:
            D_uncs[errs_too_low] = 0.03 * np.abs(Ds[errs_too_low])
        else:
            D_uncs[errs_too_low] = 0.03 * np.repeat(np.abs(Ds)[np.newaxis, :], 2, axis=0)[errs_too_low]

    return xs, Ds, D_uncs, pixel_size, window_size, pack_frac_given, pack_frac_calced