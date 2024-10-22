import common
import numpy as np
import matplotlib.pyplot as plt
import sys
import visualisation.Ds_overlapped

SHOW_TWIN_K_AXIS = False
PRESENT_SMALL = False

DISCRETE_COLORS = True

ERRORBAR_ALPHA = 0.3

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
    'f_first': 'o',
    'MSD': '_',
    'MSD_short': '_',
    'MSD_long': '_',
    'MSD_first': '_',
    'MSD_centre_of_mass_onepoint': '_',
    'MSD_centre_of_mass_proximity': '_',
    'boxcounting': '_',
    'boxcounting_first_quad': '_',
    'boxcounting_shorttime': '_',
    'boxcounting_collective': 'o',
    'timescaleint': '+',
    'timescaleint_nofit_cropped': '*',
    'D0Sk': 'o',
    'C_N_simplefit': '|',
}

def show_one_file(
        i, file, sources, PLOT_AGAINST_K, TWO_PI, logarithmic_y,
        fix_axes, filename, fig, ax,
        export_destination=None,
        show_pixel=True, show_window=True, crop_end=None
    ):

    all_Ds = []
        
    
    pack_frac_calced = None
    pack_frac_given  = None
    pixel_size = None
    window_size = None
    diameter = None

    xs = {}
    Ds = {}
    D_uncs = {}

    used_sources = []

    if DISCRETE_COLORS:
        color = ['tab:blue', 'tab:orange', 'tab:green'][i]
    else:
        color = common.colormap(i, 0, 2.5)

    # get all the data
    for source in sources:
        try:
            xs[source], Ds[source], D_uncs[source], pixel_size_temp, window_size_temp, pack_frac_given_temp, pack_frac_calced_temp, diameter_temp = visualisation.Ds_overlapped.get_L_and_D(source, file, PLOT_AGAINST_K, TWO_PI)
        
            if pixel_size_temp:       pixel_size       = pixel_size_temp
            if window_size_temp:      window_size      = window_size_temp
            if pack_frac_calced_temp: pack_frac_calced = pack_frac_calced_temp
            if pack_frac_given_temp:  pack_frac_given  = pack_frac_given_temp
            if diameter_temp:         diameter         = diameter_temp

            used_sources.append(source)

        except FileNotFoundError:
            print('FileNotFound', source)

    # if 'MSD_short' in used_sources:
    #     ax.set_ylabel('$D/D_{{MSD}}$')
    #     rescale_y = Ds['MSD_short'][0]
    # else:
    ax.set_ylabel(r'$D/D_{MSD}$')
    # ax.set_ylabel('$D$ ($\mathrm{\mu m^2/s}$)')
    #     rescale_y = 1

    D_MSD, _, _ = visualisation.Ds_overlapped.get_D0(file)
    rescale_y = D_MSD

    if diameter:
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

    assert len(xs.values()) > 0

    xmin = min([min(x) for x in xs.values() if len(x)>0]) / rescale_x
    xmax = max([max(x) for x in xs.values() if len(x)>0]) / rescale_x

    # if pack_frac_calced:
    #     pack_frac = pack_frac_calced
    # elif pack_frac_given:
    #     pack_frac = pack_frac_given
    # else:
    #     pack_frac = None
    # if pack_frac:
    #     x = (1+pack_frac)/(1-pack_frac)**3
    #     if rescale_y != 1:
    #         ax.hlines(x, xmin, xmax, label=r'$(1+\phi)/(1-\phi)^3$', color='gray')

    # do the actual plotting
    for source in used_sources:
        source_label = f'{file} {source_names.get(source, source)}'# if not source.startswith('MSD') else None
        # phi = int(file[-3:]) / 100
        # print(phi)
        # source_label = fr'$\phi={phi:.2f}$'
        # f'{file}: '+

        xs[source] /= rescale_x

        # color = colors[source]f
        ys = Ds[source] / rescale_y
        yerrs = D_uncs[source] / rescale_y

        if source in ['MSD_short', 'MSD_first']:
            ax.hlines(ys[0], xmin, xmax, color=color, linestyle='dotted', label=source_label)
            print('MSD errors hacked')
            ax.fill_between(ax.get_xlim(), ys[0]*0.97, ys[0]*1.03, facecolor=color, alpha=ERRORBAR_ALPHA)

        else:
            x_this = xs[source]
            if crop_end:
                x_this = x_this[:crop_end]
                ys     = ys    [:crop_end]
                if len(yerrs.shape) == 2:
                    yerrs = yerrs [:, :crop_end]
                else:
                    yerrs  = yerrs [:crop_end]

            ax.plot(x_this, ys, linestyle='none', marker=markers.get(source, '.'), markersize=4, color=color, label=source_label)
            ax.errorbar(x_this, ys, yerr=yerrs, linestyle='none', marker='none', alpha=ERRORBAR_ALPHA, color=color)
            # print(xs[source], ys)
        # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'
        [all_Ds.append(D) for D in ys]

            
    if not PLOT_AGAINST_K:
        if show_pixel:
            ax.vlines(pixel_size,  min(ys), max(ys), color='gray', linestyle='dotted', label='pixel size')
        if show_window:
            ax.vlines(window_size, min(ys), max(ys), color='gray', linestyle='dashed', label='window size')

    assert len(all_Ds) > 0, 'no Ds were found at all'

    if logarithmic_y:
        ax.semilogy()

    ylim_expand = 1.2
    if np.nanmax(all_Ds) - np.nanmax(all_Ds) < 0.4:
        ylim_expand = 1.5
    ymin = np.nanmin(all_Ds)/ylim_expand
    if 'MSD_short' in used_sources:
        ymin = 0.3
    ax.set_ylim(ymin, np.nanmax(all_Ds)*ylim_expand)
    if fix_axes:
        # ax.set_ylim(0.1, 5)
        pass
    # ax.set_ylabel('$D/D_0$')
    # ax.set_xticks([])
    ax.semilogx()
    if PLOT_AGAINST_K:
        

        if SHOW_TWIN_K_AXIS:
            realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
            # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
            realspace_ax.set_xlabel(r'$2\pi/k / \sigma$')
    else:
        ax.set_xlabel(r'$L / \sigma$')




def go(files, export_destination=None):

    if PRESENT_SMALL:
        figsize = (3.2, 2.8)
        figsize = (4, 3)
    else:
        figsize = (5, 4)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    for i, file in enumerate(files):

        filename = file

        show_one_file(
            i, file, [
                'boxcounting_collective',
                # 'f_first',
                # 'f_short',
                # 'f',
                # 'f_long',
                'MSD_first',
                # 'timescaleint',
                'timescaleint_nofit_cropped'
            ],
            PLOT_AGAINST_K=False, TWO_PI=True, logarithmic_y=True, fix_axes=False,
            fig=fig, ax=ax, filename=filename, show_window=False, show_pixel=False, crop_end=-4
        )

    ax.legend(fontsize=5)
    
    if export_destination:
        common.save_fig(fig, export_destination, hide_metadata=True)
    filenames = '_'.join(sys.argv[1:])
    common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_mult_{filenames}.png', dpi=200)

    
if __name__ == '__main__':
    go(sys.argv[1:])