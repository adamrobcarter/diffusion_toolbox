import countoscope_theory.nmsd
import countoscope_theory.structure_factor
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
# import sDFT_interactions
import matplotlib.cm
import visualisation.Ds_overlapped
import scipy.signal
import tqdm
import box_counting.D_of_L
import warnings
import box_counting.N_histogram

# enums
RESCALE_Y_VAR_N   = 1
RESCALE_Y_N       = 2
RESCALE_Y_PLATEAU = 3
RESCALE_Y_L2      = 4
RESCALE_X_L2 = 1
RESCALE_X_L  = 2

LINEAR_Y = False

DO_LINEAR_FIT_TO_START = False
DO_D_FROM_FIRST_POINT = True # idk why you would ever want this false

PRESENT_SMALL = False
SHOW_JUST_ONE_BOX = False

LABELS_ON_PLOT = True
LABELS_ON_PLOT_Y_SHIFT = 1.4 if SHOW_JUST_ONE_BOX else 1.25

FORCE_HIDE_LEGEND = False
SHOW_D_IN_LEGEND = False
LEGEND_LOCATION = 'upper left'

SHOW_THEORY_FIT = True
SHOW_PLATEAUS_THEORY = False
SHOW_VARIANCE = False
SHOW_MEAN = False
SHOW_PLATEAUS_OBS = False
SHOW_PLATEAU_OBS_AREA = False
SHOW_SHORT_TIME_FIT = False
SHOW_TIMESCALEINT_REPLACEMENT = False
SHOW_T_SLOPE = False

MAX_BOXES_ON_PLOT = 6
DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE = True

REVERSE_PLOT_ORDER = False
LINESTYLE = 'none'

DO_TIMESCALEINT_REPLACEMENT = False

# if SHOW_JUST_ONE_BOX:
    # LABELS_ON_PLOT = False

# RESCALE_X = RESCALE_X_L2
# RESCALE_X = RESCALE_X_L
RESCALE_X = None

# RESCALE_Y = RESCALE_Y_VAR_N
# RESCALE_Y = RESCALE_Y_PLATEAU
# RESCALE_Y = RESCALE_Y_N
RESCALE_Y = RESCALE_Y_L2
# RESCALE_Y = None

SHORTTIME_FIT_ERROR_THRESH = 0.05 # D_unc/D must be smaller than this for the point to save

TIMESCALEINT_REPLACEMENT_PLATEAU_SOURCE = 'var'

have_displayed_at_least_one = False

def go(file, ax=None, separation_in_label=False,
       linestyle='none', show_title=False,
       show_timescaleint_replacement=False, show_variance=False, labels_on_plot=True,
       rescale_x=None, rescale_y=None, legend_fontsize=7, legend_location=LEGEND_LOCATION,
       box_size_indices=None, show_nointer_theory_limits=False, max_boxes_on_plot=MAX_BOXES_ON_PLOT,
       timescaleint_replacement_plateau_source='var', nointer_theory_limit_labels=[],
       disable_ylabel=False, show_second_legend=True, show_rescaled_theory=False,
    ):
    # D0_from_fits     = [{}, {}]
    # D0_unc_from_fits = [{}, {}]
    # Dc_from_fits     = [{}, {}]
    # Dc_unc_from_fits = [{}, {}]

    if not ax:
        force_display_false = True
        # prevents anything being plotted or shown
    else:
        force_display_false = False

    LOWTIME_FIT_END = 20
    
    if rescale_x and rescale_y:
        labels_on_plot = False
        print('I disabled labels on plot as you are rescaling both axes')

    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean        = data['N2_mean']
    N2_std         = data['N2_std']
    phi            = data.get('pack_frac', np.nan)
    sigma          = data['particle_diameter']
    time_step      = data['time_step']
    depth_of_field = data.get('depth_of_field')

    box_sizes    = data['box_sizes']
    # N_mean       = data['N_mean']
    N_mean       = data.get('N_mean')
    N_var        = data.get('N_var')
    # N_var_std    = data.get('N_var_std')
    N_var_std    = np.zeros_like(N_var)
    N_var_mod    = data.get('N_var_mod')
    # N_var_losecorr = data.get('N_var_losecorr')
    N_var_losecorr    = np.zeros_like(N_var)
    # num_of_boxes = data['num_boxes']
    sep_sizes    = data['sep_sizes']

    num_timesteps = N2_mean.shape[1]
    t_all = data.get('t', np.arange(0, num_timesteps) * time_step)

    D_MSD, _, _ = visualisation.Ds_overlapped.get_D0(file)

    if REVERSE_PLOT_ORDER:
        iter = range(len(box_sizes)-1, -1, -1)
    else:
        iter = range(len(box_sizes))

    display_i = 0

    plotted_handles = []

    for box_size_index in tqdm.tqdm(iter, desc='box sizes'):
        L   = box_sizes[box_size_index]
        sep = sep_sizes[box_size_index]

        if box_size_indices:
            display = box_size_index in box_size_indices
        elif SHOW_JUST_ONE_BOX:
            display = box_size_index == 20
        elif len(box_sizes) <= max_boxes_on_plot:
            display = True
        else:
            display = box_size_index % (len(box_sizes) // max_boxes_on_plot) == 0
        
        if force_display_false:
            display = False
            # the user can provide no ax if they want to compute the Ds without plotting

        
        L_over_sigma_str = f' = {L/sigma:.2f}σ' if sigma else ''
        if display:
            print(f'L = {L:.2g}{L_over_sigma_str}')


        delta_N_sq     = N2_mean[box_size_index, :]
        delta_N_sq_err = N2_std [box_size_index, :]
        t = np.copy(t_all)
        t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()*1), 100)


        anomalous = delta_N_sq < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        assert anomalous[1] == False

        if np.any(anomalous):
            if np.sum(anomalous) > 1:
                print(f'  found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
            delta_N_sq     = delta_N_sq    [~anomalous]
            delta_N_sq_err = delta_N_sq_err[~anomalous]
            t              = t             [~anomalous]

        # assert anomalous.sum()/delta_N_sq.size < 0.8

        nmsd_nan = common.nanfrac(delta_N_sq)
        if nmsd_nan: print('nmsd nanfrac', nmsd_nan)
        nmsd_zero = np.sum(delta_N_sq==0)/delta_N_sq.size
        if nmsd_zero: print('nmsd zero', nmsd_zero)
        nmsd_negative = np.sum(delta_N_sq[1:]<0)/delta_N_sq[1:].size
        if nmsd_negative: print('nmsd negative', nmsd_negative)
        nmsd_negative = np.sum(delta_N_sq[1:]<1e-14)/delta_N_sq[1:].size
        if nmsd_negative: print('nmsd negative', nmsd_negative)
        nmsd_inf = np.sum(np.isinf(delta_N_sq))/delta_N_sq.size
        if nmsd_inf: print('nmsd inf', nmsd_inf)

        # color = matplotlib.cm.afmhot((box_size_index+2)/(len(box_sizes)+7))
        color =  common.colormap(box_size_index, 0, len(box_sizes))

        # fit to whole thing
        if True and np.isfinite(phi): # we calculate the fit even if we don't need to, because we use it for getting the angles for labels_on_plot
            if depth_of_field:
                N2_theory = lambda t, D, N: common.N2_nointer_3D(t, D, N, L, L, depth_of_field)
                type_of_fit = 'sDFT fit (no inter, 3D)'
            else:
                N_mean_for_fit = N_mean[box_size_index]
                # plateau_for_fit_mod = get_plateau(N2_mean[box_size_index, :], file)[0]

                # plateau_for_fit_mod = get_plateau(N2_mean[box_size_index, :], file)[0] * 2
                # print('aaa', get_plateau(N2_mean[box_size_index, :], file)[0], N2_mean[box_size_index, N2_mean.shape[1]//2])
                if np.isfinite(phi) and np.isfinite(sigma):
                    N2_theory = lambda t, D : countoscope_theory.nmsd.inter_2d(t, D, N_mean_for_fit, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
                    N2_theory_Lh = lambda t, D : countoscope_theory.nmsd.inter_2d_Lh(t, D, N_mean_for_fit, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma), sigma, phi)
                    type_of_fit = 'sDFT (no fit) (w/ inter.)'
                else:
                    # N2_theory = lambda t, D : countoscope_theory.nmsd.nointer_2d(t, D, N_mean_for_fit, L)
                    # type_of_fit = 'sDFT (no fit) (no inter.)'
                    pass
           
        
        
        if rescale_y:
            rescale_y_value = 1
            if rescale_y == RESCALE_Y_N:
                rescale_y_value = N_mean[box_size_index]
            elif rescale_y == RESCALE_Y_VAR_N:
                rescale_y_value = N_var[box_size_index]
            elif rescale_y == RESCALE_Y_PLATEAU:
                rescale_y_value = get_plateau(delta_N_sq, file, L, phi, sigma, D_MSD, t=t_all)[0]
            elif rescale_y == RESCALE_Y_L2:
                rescale_y_value = L**2
            delta_N_sq          /= rescale_y_value
            delta_N_sq_err      /= rescale_y_value

        if rescale_x:
            rescale_x_value = 1
            if rescale_x == RESCALE_X_L2:
                rescale_x_value = L**2
            elif rescale_x == RESCALE_X_L:
                rescale_x_value = L
            
            t_theory = t_theory / rescale_x_value
            t        = t        / rescale_x_value
        
        
        if display:
            # print('  box_size_index', box_size_index)
            # get_plateau(nmsd=N2_mean[box_size_index, :], file=file, L=L, phi=phi, sigma=sigma, display=True, t=t_all)
            
            if LINEAR_Y:
                if display_i != 0:
                    ax = ax.twinx() 
                yrange = delta_N_sq.max()
                OVERLAP = 0.5
                ylim_start =  - OVERLAP * yrange * display_i
                ylim_end = yrange * (1 + (max_boxes_on_plot - display_i)*OVERLAP)
                ax.set_ylim(ylim_start, ylim_end)
                ax.get_yaxis().set_visible(False)
            
            display_i += 1

            have_displayed_at_least_one = True


            if (rescale_x or rescale_y): # remove and False in future please
                markersize = 2
            else:
                if PRESENT_SMALL:
                    markersize = 5
                else:
                    markersize = 3

            # actual data
            if sigma and not np.isnan(sigma):
                label = rf'$L={L/sigma:.2g}\sigma$'
                if separation_in_label:
                    label += f', $s={sep/sigma:.2g}\sigma$'
            else:
                label = rf'$L={L:.2g}$'
                if separation_in_label:
                    label += f', $s={sep:.2g}$'
            if SHOW_D_IN_LEGEND:
                label += fr', $D_\mathrm{{short\:fit}}={common.format_val_and_unc(D_from_fit, D_from_fit_unc, 2)} \mathrm{{\mu m^2/s}}$'
            # ±{np.sqrt(pcov[0][0]):.3f}$'

            if DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE and delta_N_sq.size > 1000:
                # print('tt', t)
                points_to_plot = common.exponential_indices(t, 500)
                if 0 not in points_to_plot:
                    points_to_plot = np.concatenate(([0], points_to_plot))
            else:
                points_to_plot = np.index_exp[:] # this is basically a functional way of writing points_to_plot = [1:]

            data_x = t[points_to_plot]
            data_y = delta_N_sq[points_to_plot] - N2_theory(data_x, D_MSD)
            exp_plot, = ax.plot(data_x, data_y, label=label, linestyle=linestyle, marker='o', markersize=markersize, zorder=-1, color=color)
            # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:]/np.sqrt(num_of_boxes[box_size_index]), label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            plotted_handles.append(exp_plot)


    if not force_display_false:
        assert have_displayed_at_least_one, 'display was false for all L'

        legend = ax.legend(handles=plotted_handles, fontsize=legend_fontsize, loc=legend_location)
        # common.set_legend_handle_size(legend)
        
        ax.semilogx()
        if rescale_x == RESCALE_X_L:
            xlabel = '$t/L ($\mathrm{s/\mu m})$'
        elif rescale_x == RESCALE_X_L2:
            xlabel = '$t/L^2$ ($\mathrm{s/\mu m^2}$)'
        else:
            xlabel = '$t$ ($\mathrm{s}$)'
        if rescale_y == None:
            ylabel = r'$\langle \Delta N^2(t) \rangle$ ($\mathrm{\mu m^2}$)'
        elif rescale_y == RESCALE_Y_N:
            ylabel = r'$\langle \Delta N^2(t) \rangle / \langle N\rangle$'
        elif rescale_y == RESCALE_Y_VAR_N:
            ylabel = r'$\langle \Delta N^2(t) \rangle / Var(N)$'
        elif rescale_y == RESCALE_Y_PLATEAU:
            ylabel = r'$\langle \Delta N^2(t) \rangle / \mathrm{plateau}$'
        elif rescale_y == RESCALE_Y_L2:
            ylabel = r'$\langle \Delta N^2(t) \rangle / L^2$ ($\mathrm{\mu m^{-2}}$)'
        ax.set_xlabel(xlabel)
        if not disable_ylabel: ax.set_ylabel(ylabel)

        title = file
        # title = f'Simulated colloids in RCP spheres\n$\phi={phi:.3f}$'
        if not np.isnan(phi):
            title += f', $\phi_\mathrm{{calc}}={phi:.3f}$'
        if not np.isnan(sigma):
            title += f', $\sigma={sigma:.2f}\mathrm{{\mu m}}$'
        if sigma_calced := data.get('particle_diameter_calced'):
            title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
            # print('sigma calced hidden from legend')
        if show_title:
            ax.set_title(title)


if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'counted_'):
        
        figsize = (6, 4.5)
        if PRESENT_SMALL:
            figsize = (4.5, 4)
            figsize = (3.5, 3.2)
            
        fig, ax = plt.subplots(1, 1, figsize=figsize)
            
        go(file,
           ax=ax,
           linestyle=LINESTYLE,
           show_variance=SHOW_VARIANCE,
           show_timescaleint_replacement=SHOW_TIMESCALEINT_REPLACEMENT,
           timescaleint_replacement_plateau_source=TIMESCALEINT_REPLACEMENT_PLATEAU_SOURCE,
           show_title=not PRESENT_SMALL,
           labels_on_plot=LABELS_ON_PLOT,
           rescale_x=RESCALE_X, rescale_y=RESCALE_Y,
           max_boxes_on_plot=MAX_BOXES_ON_PLOT,
           legend_fontsize=10
        )
        
        common.save_fig(fig, f'box_counting/figures_png/nmsd_minus_theory_{file}.png', dpi=200)