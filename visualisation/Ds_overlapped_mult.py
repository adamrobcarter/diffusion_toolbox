import common
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys
import visualisation.Ds_overlapped

SHOW_TWIN_K_AXIS = False
PRESENT_SMALL = False
PLOT_AGAINST_K = False

DISCRETE_COLORS = True

ERRORBAR_ALPHA = 0.3
LEGEND_FONTSIZE = 6


SOURCES = [
    'f_first_first',
    # 'D0Sk_theory',
    # 'f_t0.5',
    # 'f_t2',
    # 'f_t8',
    # 'f_t32',
    # 'f_t256',
    # 'f_t1024',
    
    # 'f_t128',
    # 'f_t512',
    # 'f_t2048',

    # 'F_first32_first',
    # 'f_first',
    # 'f_short',
    # 'f',
    # 'f_long',
    # 'MSD_first',
    # 'boxcounting_collective_var',
    # 'timescaleint_var',
    # 'timescaleint_fixexponent_var',
    # 'timescaleint_nmsdfitinter'

    'timescaleint_nofit_cropped_var',

    # 'D_of_L_theory',
    # 'D0Sk_theory'
]

source_names = {
    'DDM': 'DDM',
    'DDM_short': 'DDM short',
    'DDM_long': 'DDM long',
    'f': '$f(k, t)$',
    'Fs': '$F_s(k, t)$',
    'f_short': '$f(k, t)$ short',
    'Fs_short': '$F_s(k, \mathrm{short})$',
    'f_long': '$f(k, t)$ long',
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

marker_index = {
    'DDM_short': '*',
    'f_short': 'x',
    'DDM_long': '*',
    'f_long': '_',
    'DDM': '*',
    'f': '*',
    'f_first': 'x',
    'f_first_first': 'x',
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
    'timescaleint_var': '+',
    'timescaleint_nofit_cropped_var': '*',
    'D0Sk': 'o',
    'C_N_simplefit': '|',
    'D_of_L_theory': 'none',
    'D0Sk_theory': 'none',
}

linestyles = {
    'D0Sk_theory': '-',
    'D_of_L_theory': '-',
    'D_of_L_theory_Lh': '-',
    'dominiguez_theory': '-',
}

def show_one_file_and_source(
        file, source, PLOT_AGAINST_K, TWO_PI, logarithmic_y,
        ax,
        show_pixel=True, show_window=True, crop_end=None,
        allow_rescale_y=True, linestyle=None,
        label=None, errorbar_alpha=ERRORBAR_ALPHA,
        marker=None, color=None, disable_ylabel=False,
        fade_out_thresh=np.inf, fade_out_alpha=0.5,
    ):
    
    try:
        D_MSD, sigma, phi = visualisation.Ds_overlapped.get_D0(file)
    except:
        print('RECONSIDER THIS')
        D_MSD = 1
        sigma = None
        phi = None
        allow_rescale_y = False

    diameter = sigma

    # get all the data
    # for source in sources:
    xs, Ds, D_uncs, pixel_size, window_size, pack_frac_given, pack_frac_calced = visualisation.Ds_overlapped.get_L_and_D(source, file, PLOT_AGAINST_K, TWO_PI, D_MSD=D_MSD, sigma=sigma, phi=phi)


    if allow_rescale_y:
        try:
            rescale_y = D_MSD
            if not disable_ylabel: ax.set_ylabel(r'$D/D_\mathrm{self}$')
            no_rescale = False
        except FileNotFoundError as err:
            no_rescale = True
            print(err)
    else:
        no_rescale = True
    
    if no_rescale:
        rescale_y = 1
        if not disable_ylabel: ax.set_ylabel(r'$D$ ($\mathrm{\mu m^2/s}$)')

    if diameter:
        if PLOT_AGAINST_K:
            rescale_x = 1/diameter
            ax.set_xlabel(r'$k \sigma$')
        else:
            rescale_x = diameter
            ax.set_xlabel(r'$L/\sigma$')
    else:
        # assert False
        rescale_x = 1
        if PLOT_AGAINST_K:
            ax.set_xlabel(r'$k$')
        else:
            ax.set_xlabel(r'$L$')

    assert len(xs) > 0, f'xs was empty for {file}'
    xmin = min(xs) / rescale_x
    xmax = max(xs) / rescale_x
    # print('xlim', xmin, xmax)
    # ax.set_xlim(xmin, xmax)

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
    if not label:
        label = f'{file} {source_names.get(source, source)}'# if not source.startswith('MSD') else None

    xs /= rescale_x
    ys = Ds / rescale_y
    yerrs = D_uncs / rescale_y

    if source in ['MSD_short', 'MSD_first']:
        ax.hlines(ys[0], xmin, xmax, color=color, linestyle='dotted', label=label)
        print('MSD errors hacked')
        ax.fill_between(ax.get_xlim(), ys[0]*0.97, ys[0]*1.03, facecolor=color, alpha=errorbar_alpha)

    else:
        if crop_end:
            print('cropping', crop_end)
            xs = xs[:crop_end]
            ys     = ys    [:crop_end]
            if len(yerrs.shape) == 2:
                yerrs = yerrs [:, :crop_end]
            else:
                yerrs  = yerrs [:crop_end]
        zorder = -1 if 'theory' in source else 0
        if not linestyle:
            linestyle = linestyles.get(source, 'none') # try 'None' if this errors out
            
        if not marker:
            marker = marker_index.get(source, 'o')
        if 'theory' in source:
            marker = 'None'
        
        # print(source, 'linestyles', use_linestyle, 'marker', marker)
        # print(file, source, 'plotting', ys.size, ys)
        fade_out_thresh = np.inf
        not_faded = xs < fade_out_thresh
        # THIS NEEDS TO BE PER SOURCE I THINK, A FULL 2D ARRAY
        faded     = xs > fade_out_thresh
        ax.plot(xs[not_faded], ys[not_faded], linestyle=linestyle, marker=marker, markersize=4, color=color, label=label, zorder=zorder, linewidth=1)
        ax.plot(xs[faded    ], ys[faded    ], linestyle=linestyle, marker=marker, markersize=4, color=color,                          zorder=zorder, linewidth=1, alpha=fade_out_alpha)
        ax.errorbar(xs, ys, yerr=yerrs, linestyle='none', marker='None', alpha=errorbar_alpha, color=color, zorder=zorder)
        # print(xs[source], ys)

        log_y = np.log10(ys)
        err = np.sum((log_y[1:]-log_y[:-1])**2)
        # print('sum err', file, err)

    # assert not np.any(np.isnan(Ds)), 'nan was found in Ds'=

            
    if not PLOT_AGAINST_K:
        if show_pixel:
            ax.vlines(pixel_size,  min(ys), max(ys), color='gray', linestyle='dotted', label='pixel size')
        if show_window:
            ax.vlines(window_size, min(ys), max(ys), color='gray', linestyle='dashed', label='window size')

    assert len(Ds) > 0, 'no Ds were found at all'
    assert np.isfinite(Ds).any(), 'Ds were found but they were all nan'

    if logarithmic_y:
        ax.semilogy()
        ymin, ymax = ax.get_ylim()
        if ymax/ymin < 50:
            ax.yaxis.set_minor_formatter(matplotlib.ticker.LogFormatter()) # prevent scientific notation on axes
            ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatter()) # prevent scientific notation on axes

    ylim_expand = 1.5
    ymin = max(0.01, np.nanmin(ys[ys > 0])/ylim_expand)
    print('ymin', ymin, np.nanmin(ys[ys > 0]))
    ymax = np.nanquantile(ys, 0.95)*ylim_expand*1.2
    if source == 'MSD_short':
        ymin = 0.3
    ax.set_ylim(ymin, ymax)
    
    ax.semilogx()
    if PLOT_AGAINST_K:

        if SHOW_TWIN_K_AXIS:
            realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
            # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
            realspace_ax.set_xlabel(r'$2\pi/k / \sigma$')
    else:
        ax.set_xlabel(r'$L / \sigma$')


def go(files_and_sources, ax, plot_against_k=False, legend_fontsize=None,
       discrete_colors=False, logarithmic_y=False, labels=None,
       errorbar_alpha=ERRORBAR_ALPHA, markers=None, source_labels=None,
       allow_rescale_y=True, colors=None, linestyles=None, disable_ylabel=False,
       fade_out_alpha=0.5, show_Dcoll=False, allow_missing_files=False):
    # colors can be len(files) or len(files) x len(sources)
    
    if colors:
        assert len(colors) == len(files_and_sources)
    if labels:
        assert len(labels) == len(files_and_sources)
    if linestyles:
        assert len(linestyles) == len(files_and_sources)

    for i, (file, source) in enumerate(files_and_sources):
        if labels:
            label = labels[i]
        else:
            label = None

        if type(markers) == str:
            marker = markers
        elif type(markers) == list:
            marker = markers[i]
        else:
            marker = None

        if linestyles:
            linestyle = linestyles[i]
        else:
            linestyle = None

        if colors:
            color = colors[i]
        else:
            if discrete_colors:
                color = common.tab_color(i)
            else:
                color = common.colormap(i, 0, len(files_and_sources))
        
        try:
            show_one_file_and_source(
                file, source,
                PLOT_AGAINST_K=plot_against_k, TWO_PI=True, logarithmic_y=logarithmic_y,
                ax=ax, show_window=False, show_pixel=False,
                allow_rescale_y=allow_rescale_y, linestyle=linestyle,
                label=label, errorbar_alpha=errorbar_alpha,
                marker=marker, color=color,
                disable_ylabel=disable_ylabel,
                fade_out_thresh=None, fade_out_alpha=fade_out_alpha,
            )
            
        except FileNotFoundError as err:
            if allow_missing_files:
                print('FileNotFound', err)
            else:
                raise err

    ax.hlines(1, *ax.get_xlim(), linestyle=(0, (0.7, 0.7)), color='darkgray')

    legend_margin = -0.015
    ax.legend(
        fontsize=legend_fontsize,
        loc='upper left' if not plot_against_k else 'upper right',
        bbox_to_anchor=(legend_margin, legend_margin, 1-2*legend_margin, 1-2*legend_margin)
    )
    
if __name__ == '__main__':
    
    if PRESENT_SMALL:
        figsize = (3.2, 2.8)
        figsize = (4, 3)
    else:
        figsize = (5, 4)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    files = common.files_from_argv('isf/data/', 'F_first_')


    files_and_sources = []
    [[files_and_sources.append((file, source)) for source in SOURCES] for file in files]
    colors = []
    [[colors.append(common.tab_color(i)) for source in SOURCES] for i in range(len(files))]
    go(
        files_and_sources,
        ax,
        linestyles=['-']*len(files_and_sources),
        legend_fontsize=LEGEND_FONTSIZE,
        # discrete_colors=DISCRETE_COLORS,
        colors=colors,
        allow_rescale_y=True,
        plot_against_k=PLOT_AGAINST_K,
        allow_missing_files=True
    )
    # ax.set_ylim(0.9, 2.5)
    
    filenames = '_'.join(files)
    common.save_fig(fig, f'visualisation/figures_png/Ds_overlapped_mult_{filenames}.png', dpi=200)
