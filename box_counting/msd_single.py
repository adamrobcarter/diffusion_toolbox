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


# enums
RESCALE_Y_VAR_N   = 1
RESCALE_Y_N       = 2
RESCALE_Y_PLATEAU = 3
RESCALE_X_L2 = 1
RESCALE_X_L  = 2

LINEAR_Y = False

DO_LINEAR_FIT_TO_START = False
DO_D_FROM_FIRST_POINT = True # idk why you would ever want this false

PRESENT_SMALL = False
SHOW_JUST_ONE_BOX = False

LABELS_ON_PLOT = False
LABELS_ON_PLOT_Y_SHIFT = 1.4 if SHOW_JUST_ONE_BOX else 1.25

FORCE_HIDE_LEGEND = False
SHOW_D_IN_LEGEND = False
LEGEND_LOCATION = 'upper left'

SHOW_THEORY_FIT = False
SHOW_PLATEAUS_THEORY = False
SHOW_VARIANCE = False
SHOW_MEAN = False
SHOW_PLATEAUS_OBS = True
SHOW_PLATEAU_OBS_AREA = False
SHOW_SHORT_TIME_FIT = False
SHOW_TIMESCALEINT_REPLACEMENT = True
SHOW_T_SLOPE = False

MAX_BOXES_ON_PLOT = 6
DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE = False

REVERSE_PLOT_ORDER = False
LINESTYLE = '-'

# if SHOW_JUST_ONE_BOX:
    # LABELS_ON_PLOT = False

# RESCALE_X = RESCALE_X_L2
# RESCALE_X = RESCALE_X_L
RESCALE_X = None

# RESCALE_Y = RESCALE_Y_VAR_N
# RESCALE_Y = RESCALE_Y_PLATEAU
# RESCALE_Y = RESCALE_Y_N
RESCALE_Y = False

if RESCALE_X and RESCALE_Y:
    LABELS_ON_PLOT = False

SHORTTIME_FIT_ERROR_THRESH = 0.05 # D_unc/D must be smaller than this for the point to save

figsize = (6, 4.5)
if PRESENT_SMALL:
    figsize = (4.5, 4)
    figsize = (3.5, 3.2)

have_displayed_at_least_one = False

def get_plateau_range(file):
    if 'eleanorlong' in file: # hacky, pls don't do this
        start_index = -70000
        end_index   = -20000
        if file.startswith('eleanorlong066'):
            start_index = -100000
        if file.startswith('eleanorlong001'):
            start_index = -100000
            end_index   = -70000
    elif 'brennan' in file: # hacky, pls don't do this
        if '034' in file:
            start_index = -20000
            end_index   = -10000
        elif '066' in file:
            start_index = -10000
            end_index   = - 5000
    elif 'marine3' in file: # hacky, pls don't do this
        start_index = -1900
        end_index   = -1000
    else: # used to be -300, -100
        start_index = -600
        end_index   = -400

    return start_index, end_index

gradfig, (gradax1, gradax2, gradax3, gradax4) = plt.subplots(4, 1, figsize=(5, 4*4))
gradax1.semilogx()
gradax2.loglog()
gradax3.loglog()
gradax4.loglog()

def get_plateau(nmsd, file, L, phi, sigma, t, display=False):

    if display:
        assert nmsd.size == t.size, f'nmsd {nmsd.shape}, t {t.shape}'
        grad = np.gradient(nmsd, t)
        gradax1.scatter(t, grad, s=1)
        gradax2.scatter(t, np.abs(grad), s=1)

        t_indexes = np.arange(0, t.size, 1)

        width = 0.05
        filter_size = [int(s*width) for s in t_indexes]
        # print(filter_size[:100])
        # print(filter_size[::500])
        # nmsd2 = [np.abs(grad)[i-filter_size[i]:min(i+filter_size[i]+1, t.size-1)].mean() for i in t_indexes]
        # gradax3.scatter(t, nmsd2, s=1)

        points = common.exponential_integers(1, t.size-1, 100)
        # gradax4.scatter(t[points], np.abs(grad)[points], s=1)
        nmsd3 = [nmsd[points[i]:points[i+1]].mean() for i in range(points.size-1)]
        # if you change grad to MSD it looks clean AF...
        gradax4.scatter(t[points][:-1], nmsd3, s=1)

        # gradax3.scatter(t, scipy.signal.savgol_filter(np.abs(grad)[1000:], 1000, 3))
        

    start_index, end_index = get_plateau_range(file)  

    used_data = nmsd[start_index:end_index] # used to be -300:-100, we could do with a more inteligent method (use the gradient (smoothed?))
    # used_data = nmsd[-300:-100]
    assert used_data.size > 10, f'used_data.shape == {used_data.shape}'

    assert end_index-start_index > 1000
    num_points_to_mean_over = int((end_index-start_index)/10)
    # print(num_points_to_mean_over)

    points_to_mean_over = common.exponential_integers(1, used_data.size-1, num_points_to_mean_over)
    # print('points', points_to_mean_over)

    obs_plat_mean = used_data[points_to_mean_over].mean()
    obs_plat_std = used_data[points_to_mean_over].std()
    

    if file.startswith('eleanorlong001') or file == 'eleanor0.01':
        D0 = 0.042
    elif file.startswith('eleanorlong034'):
        D0 = 0.031
    else:
        raise Exception()

    if False:
        nmsd_th = countoscope_theory.nmsd.inter_2d(t, D0, 1, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
        nmsd_th_plat = countoscope_theory.nmsd.plateau_inter_2d(1, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
        nmsd_th /= nmsd_th_plat
        print('  last nmsd', nmsd_th[-10:])
        points = common.exponential_integers(1, t.size-1, t.size//10)
        func = lambda plat : np.sum((plat*nmsd_th[points] - nmsd[points])**2)
        # popt, pcov = scipy.optimize.curve_fit(func, [int(i) for i in np.arange(0, t.size)], nmsd)
        res = scipy.optimize.minimize(func, x0=[obs_plat_mean])
        assert res.success
        plat_optimised = res.x[0]
        # print('got', popt)
    
        print('  plats', obs_plat_mean, plat_optimised, obs_plat_mean/plat_optimised)

    return obs_plat_mean, obs_plat_std

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'counted_'):

        # D0_from_fits     = [{}, {}]
        # D0_unc_from_fits = [{}, {}]
        # Dc_from_fits     = [{}, {}]
        # Dc_unc_from_fits = [{}, {}]

        LOWTIME_FIT_END = 20

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

        data = common.load(f'box_counting/data/counted_{file}.npz')
        # data = common.load(f'data/counted_driftremoved_{phi}.npz')
        N2_mean        = data['N2_mean']
        N2_std         = data['N2_std']
        phi            = data['pack_frac']
        sigma          = data['particle_diameter']
        time_step      = data['time_step']
        depth_of_field = data.get('depth_of_field')

        # N_stats        = data['N_stats']
        # box_sizes    = N_stats[:, 0]
        # N_mean       = N_stats[:, 1]
        # N_var        = N_stats[:, 2]
        # num_of_boxes = N_stats[:, 5]
        box_sizes    = data['box_sizes']
        # N_mean       = data['N_mean']
        N_mean       = data.get('N_mean')
        N_var        = data.get('N_var')
        N_var_mod    = data.get('N_var_mod')
        # num_of_boxes = data['num_boxes']
        sep_sizes    = data['sep_sizes']

        print('L/sigma', box_sizes/sigma)

        num_timesteps = N2_mean.shape[1]
        num_boxes     = N2_mean.shape[0]
        t_all = np.arange(0, num_timesteps) * time_step

        start_index, end_index = get_plateau_range(file)
        if not LINEAR_Y:
            ax.vlines(t_all[start_index], 1e-3, 1e3, alpha=0.5)
            ax.vlines(t_all[end_index], 1e-3, 1e3, alpha=0.5)

        # reduce = 1
        # t        = t      [::reduce]
        # # t_theory = np.logspace()
        # N2_mean  = N2_mean[:, ::reduce]
        # N2_std   = N2_std [:, ::reduce]
        

        Ds_for_saving = []
        D_uncs_for_saving = []
        Ls_for_saving = []

        Ds_shorttime_for_saving = []
        D_uncs_shorttime_for_saving = []
        Ls_shorttime_for_saving = []

        Ds_first_quad_for_saving = []
        D_uncs_first_quad_for_saving = []
        Ls_first_quad_for_saving = []
        
        Ds_for_saving_collective = []
        D_uncs_for_saving_collective = []
        Ls_for_saving_collective = []

        if REVERSE_PLOT_ORDER:
            iter = range(len(box_sizes)-1, -1, -1)
        else:
            iter = range(len(box_sizes))

        display_i = 0

        for box_size_index in iter:
            L   = box_sizes[box_size_index]
            sep = sep_sizes[box_size_index]
            
            L_over_sigma_str = f' ={L/sigma:.2f}σ' if sigma else ''
            print(f'L={L:.2g}{L_over_sigma_str} ({box_size_index}/{len(box_sizes)})')


            delta_N_sq     = N2_mean[box_size_index, :]
            delta_N_sq_err = N2_std [box_size_index, :]
            t = np.copy(t_all)
            t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()*1), 100)

            anomalous = delta_N_sq < 1e-14
            anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
            if np.any(anomalous):
                print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
                delta_N_sq     = delta_N_sq    [~anomalous]
                delta_N_sq_err = delta_N_sq_err[~anomalous]
                t              = t             [~anomalous]
            assert anomalous.sum()/delta_N_sq.size < 0.8
            
            # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
            # N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - common.famous_f(4*D0*t/L**2) * common.famous_f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

            # fit_func = N2_func_full
            # popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
            # D0 = popt[0]
            # r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))


            # color = matplotlib.cm.afmhot((box_size_index+2)/(len(box_sizes)+7))
            color =  common.colormap(box_size_index, 0, len(box_sizes))

            # computed theory interactions
            # D0 = { # countoscope paper, table 1
            #     'alice0.02': 0.0416,
            #     'alice0.02_overlapped': 0.0416,
            #     'alice0.34': 0.0310,
            #     'alice0.66': 0.0175
            # }[file]

            # fit to whole thing
            if SHOW_THEORY_FIT:
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
                log_N2_theory = lambda t, *args : np.log(N2_theory(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
                log_N2_theory_Lh = lambda t, *args : np.log(N2_theory_Lh(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
                
                fitting_points = common.exponential_integers(1, t.size//2)
                # p0 = (0.05, N_mean[box_size_index])
                p0 = [0.05]
                popt, pcov = scipy.optimize.curve_fit(log_N2_theory, t[fitting_points], np.log(delta_N_sq[fitting_points]), p0=p0, maxfev=2000)
                popt_Lh, pcov_Lh = scipy.optimize.curve_fit(log_N2_theory_Lh, t[fitting_points], np.log(delta_N_sq[fitting_points]), p0=p0, maxfev=2000)
                # D_from_fit = popt[0]
                # D_from_fit_unc = np.sqrt(pcov[0][0])

                # D0, _, _ = visualisation.Ds_overlapped.get_D0(file)
                # popt = (D0,)
                # print('D0 was', D0)
                
                N2_theory_points = N2_theory(t_theory, *popt)
                N2_theory_points_Lh = N2_theory_Lh(t_theory, *popt_Lh)
                print('  plats', get_plateau(N2_mean[box_size_index, :], file, L, phi, sigma, t=t_all)[0], N2_theory_points[-1])

            
            # fit to whole thing 2 - replace timescaleint
            def timescaleint_replacement(nmsd):
                plateau, plateau_unc = get_plateau(nmsd, file, L, phi, sigma, t=t_all)
                print(f'  plat unc {plateau_unc/plateau:.3f}')
                N2_theory2 = lambda t, D : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
                log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
                
                fitting_points2 = common.exponential_integers(1, t.size//2)
                # p0 = (0.05, N_mean[box_size_index])
                p02 = [0.05]
                popt2, pcov2 = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points2], np.log(nmsd[fitting_points2]), p0=p02, maxfev=2000)
                # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
                D_from_fit2 = popt2[0]
                D_from_fit_unc2 = np.sqrt(pcov2[0][0])

                if D_from_fit_unc2/D_from_fit2 < plateau_unc/plateau:
                    D_from_fit_unc2 = D_from_fit2 * plateau_unc/plateau
                    # print(f'  increased error to {plateau_unc/plateau:.3f}')

                fit_ys = N2_theory2(t_theory, *popt2)

                return D_from_fit2, D_from_fit_unc2, fit_ys
            def timescaleint_replacement_fitplateau(nmsd):
                N2_theory2 = lambda t, D, plateau : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
                log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
                
                fitting_points2 = common.exponential_integers(1, t.size//2)
                # p0 = (0.05, N_mean[box_size_index])
                p02 = [0.05, 2*N_var[box_size_index]]
                popt2, pcov2 = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points2], np.log(nmsd[fitting_points2]), p0=p02, maxfev=2000)
                # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
                D_from_fit2 = popt2[0]
                D_from_fit_unc2 = np.sqrt(pcov2[0][0])

                fit_ys = N2_theory2(t_theory, *popt2)

                return D_from_fit2, D_from_fit_unc2, fit_ys

            tsi_replace_func = timescaleint_replacement
            # tsi_replace_func = timescaleint_replacement_fitplateau
            Dc, Dc_unc, tsi_replacement_ys = tsi_replace_func(N2_mean[box_size_index, :])
            print(f'  timescaleint replacement D={common.format_val_and_unc(Dc, Dc_unc, latex=False)}')
            Dc_lower, Dc_unc_lower, _ = tsi_replace_func(N2_mean[box_size_index, :]-N2_std[box_size_index])
            Dc_upper, Dc_unc_upper, _ = tsi_replace_func(N2_mean[box_size_index, :]+N2_std[box_size_index])
            # Dc_unc_final = max(Dc_unc, abs(Dc-Dc_lower), abs(Dc-Dc_upper))
            Dc_unc_final = Dc_unc # errors are way overestimated with the above line
            print(f'  final Dc unc {Dc_unc_final/Dc:.3f}')

            Ds_for_saving_collective.append(Dc)
            D_uncs_for_saving_collective.append(Dc_unc_final)
            Ls_for_saving_collective.append(L)
            
            if RESCALE_Y:
                if RESCALE_Y == RESCALE_Y_N:
                    rescale = N_mean[box_size_index]
                elif RESCALE_Y == RESCALE_Y_VAR_N:
                    rescale = N_var[box_size_index]
                elif RESCALE_Y == RESCALE_Y_PLATEAU:
                    rescale = get_plateau(delta_N_sq, file, L, phi, sigma, t=t_all)[0]
                delta_N_sq       /= rescale
                delta_N_sq_err   /= rescale
                if SHOW_THEORY_FIT:
                    N2_theory_points    /= rescale
                    N2_theory_points_Lh /= rescale

            if RESCALE_X:
                if RESCALE_X == RESCALE_X_L2:
                    t /= L**2
                elif RESCALE_X == RESCALE_X_L:
                    rescale = L**0.3
                    
                t_theory /= rescale
                t        /= rescale
            
            info = fr'L = {L/sigma:.1f}σ,'
            
            # linear fit to start
            if DO_LINEAR_FIT_TO_START:
                fit_end = 4
                fit_func_2 = lambda t, D : 8 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(t * D / L**2)
                
                popt, pcov = scipy.optimize.curve_fit(fit_func_2, t[1:fit_end], delta_N_sq[1:fit_end])
                D_from_shorttime = popt[0]
                D_unc_from_shorttime = np.sqrt(pcov)[0, 0]
                shorttime_fit_is_good = D_unc_from_shorttime/D_from_shorttime < SHORTTIME_FIT_ERROR_THRESH
                
                if shorttime_fit_is_good:
                    print(f'  short good D_unc/D={D_unc_from_shorttime/D_from_shorttime:.3f}')
                    Ds_shorttime_for_saving.append(D_from_shorttime)
                    D_uncs_shorttime_for_saving.append(D_unc_from_shorttime)
                    Ls_shorttime_for_saving.append(L)
                else:
                    print(f'  short bad D_unc/D={D_unc_from_shorttime/D_from_shorttime:.3f}')

            
            # quadratic fit to start
            # fit_end = 6
            # fit_func_quad = lambda t, D : 8 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(t * D / L**2)
            
            # popt, pcov = scipy.optimize.curve_fit(fit_func_quad, t[1:fit_end], delta_N_sq[1:fit_end])
            # D_from_shorttime_quad = popt[0]
            # D_unc_from_shorttime_quad = np.sqrt(pcov)[0, 0]
            
            # if D_unc_from_shorttime_quad/D_from_shorttime_quad > 0.03:
            #     pass
            # else:
            if DO_D_FROM_FIRST_POINT:
                D_from_first_quad = np.pi * L**2 / ( 4 * time_step ) * (1 - np.sqrt(1 - delta_N_sq[1]/(2*N_mean[box_size_index])))**2
                # error propagation for that is complicated (type d/dx A(1-sqrt(1-x/(2B)))^2 into wolfram so let's hack)
                D_unc_from_first_quad = delta_N_sq_err[1] / delta_N_sq[1] * D_from_first_quad
                Ds_first_quad_for_saving.append(D_from_first_quad)
                D_uncs_first_quad_for_saving.append(D_unc_from_first_quad)
                Ls_first_quad_for_saving.append(L)
            

            if len(box_sizes) <= MAX_BOXES_ON_PLOT:
                display = True
            else:
                display = box_size_index % (len(box_sizes) // MAX_BOXES_ON_PLOT) == 0

            if SHOW_JUST_ONE_BOX:
                display = box_size_index == 20

            if display:
                get_plateau(N2_mean[box_size_index, :], file, L, phi, sigma, display=True, t=t_all)
                
                if LINEAR_Y:
                    if display_i != 0:
                        ax = ax.twinx() 
                    yrange = delta_N_sq.max()
                    ylim_range = yrange * MAX_BOXES_ON_PLOT
                    OVERLAP = 0.5
                    ylim_start =  - OVERLAP * yrange * display_i
                    ylim_end = yrange * (1 + (MAX_BOXES_ON_PLOT - display_i)*OVERLAP)
                    ax.set_ylim(ylim_start, ylim_end)
                    ax.get_yaxis().set_visible(False)
                
                display_i += 1

                have_displayed_at_least_one = True

                if SHOW_MEAN:
                    ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color=color, linewidth=1, linestyle='dashdot', label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
                if SHOW_VARIANCE:
                    ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$2\mathrm{Var}(N)$' if box_size_index==0 else None)
                    ax.hlines(2*N_var_mod[box_size_index], t.min(), t.max(), linestyles='dotted', color='grey', linewidth=1, label=r'$2\mathrm{Var}(N)$' if box_size_index==0 else None)


                if not RESCALE_Y and SHOW_THEORY_FIT:
                    ax.plot(t_theory[1:], N2_theory_points[1:], color='gray', linewidth=1, label=type_of_fit if box_size_index==0 else None)
                    ax.plot(t_theory[1:], N2_theory_points_Lh[1:], color='white', linewidth=1, label=type_of_fit if box_size_index==0 else None)
                # label += fr', $D_\mathrm{{fit}}={popt[0]:.3g}\pm {np.sqrt(pcov[0][0]):.3g} \mathrm{{\mu m^2/s}}$'

                if (RESCALE_X or RESCALE_Y): # remove and False in future please
                    markersize = 2
                else:
                    if PRESENT_SMALL:
                        markersize = 5
                    else:
                        markersize = 3

                # actual data
                if LABELS_ON_PLOT:
                    label = label='observations' if box_size_index==0 else ''
                else:
                    if sigma and not np.isnan(sigma):
                        label = rf'$L={L/sigma:.2g}\sigma$'
                        label += f', $s={sep/sigma:.2g}\sigma$'
                    else:
                        label = rf'$L={L:.2g}$'
                        label += f', $s={sep:.2g}$'
                if SHOW_D_IN_LEGEND:
                    label += fr', $D_\mathrm{{short\:fit}}={common.format_val_and_unc(D_from_fit, D_from_fit_unc, 2)} \mathrm{{\mu m^2/s}}$'
                # ±{np.sqrt(pcov[0][0]):.3f}$'
                print('  ', delta_N_sq.size, common.nanfrac(delta_N_sq))

                if DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE and delta_N_sq.size > 1000:
                    points_to_plot = common.exponential_integers(1, delta_N_sq.size-1, 500)
                else:
                    points_to_plot = np.index_exp[1:] # this is basically a functional way of writing points_to_plot = [1:]
                
                exp_plot = ax.plot(t[points_to_plot], delta_N_sq[points_to_plot], label=label, linestyle=LINESTYLE, marker='o', markersize=markersize, zorder=-1, color=color)
                # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:]/np.sqrt(num_of_boxes[box_size_index]), label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
                # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            
                if LABELS_ON_PLOT:
                    t_index_for_text = int(t_theory.size // 1.6)
                    angle = np.tan(np.gradient(N2_theory_points, t_theory)[t_index_for_text]) * 180/np.pi
                    # plt.scatter(t_theory[t_index_for_text], N2_theory_points[t_index_for_text])
                    
                    if sigma and not np.isnan(sigma):
                        L_label = rf'$L={L/sigma:.2g}\sigma$'
                    else:
                        L_label = rf'$L={L:.2g}$'
                    if SHOW_JUST_ONE_BOX:
                        ax.text(t_theory[t_index_for_text+6], N2_theory_points[t_index_for_text+6]/LABELS_ON_PLOT_Y_SHIFT, L_label,
                                horizontalalignment='center', color=color, fontsize=9)
                    else:
                        ax.text(t_theory[t_index_for_text+6], N2_theory_points[t_index_for_text+6]*LABELS_ON_PLOT_Y_SHIFT, L_label,
                                horizontalalignment='center', color=color, fontsize=9,
                                transform_rotates_text=True, rotation=angle, rotation_mode='anchor')
                
                # linear fit to start
                if DO_LINEAR_FIT_TO_START:
                    if not shorttime_fit_is_good: # was 0.03
                        print(f'  skipping short time fit at L={L}um, D_unc/D={D_unc_from_shorttime/D_from_shorttime:.2f}')
                    else:
                        print(f'  short time fit: D={D_from_shorttime:.4f}')
                        # D_ratio = D_from_shorttime/D_from_fit
                        # # print(f'D_short / D_fit = {D_ratio:.2f}')
                        # if D_ratio > 1.5 or 1/D_ratio > 1.5:
                        #     print(f'problem! D fit = {common.format_val_and_unc(D_from_fit, D_from_fit_unc, 2)} D shorttime = {common.format_val_and_unc(D_from_shorttime, D_unc_from_shorttime, 2)}')
                        if SHOW_SHORT_TIME_FIT:
                            ax.plot(t[1:fit_end], fit_func_2(t[1:fit_end], *popt), linestyle=':', color=common.FIT_COLOR, linewidth=2)

                if DO_D_FROM_FIRST_POINT:
                    print(f'  from first point quad: D={D_from_first_quad:.4f}')

                if SHOW_T_SLOPE:
                    t_half_scaling_line_offset = 5
                    x1, y1 = t[1]*t_half_scaling_line_offset, delta_N_sq[1]
                    t_half_scaling_line_size = 5
                    x2, y2 = x1*t_half_scaling_line_size, y1*np.sqrt(t_half_scaling_line_size)
                    ax.plot([x1, x2], [y1, y2], color=common.FIT_COLOR)
                    ax.text(x2, y1, '$t^{1/2}$', ha='right', color=common.FIT_COLOR)

                if SHOW_TIMESCALEINT_REPLACEMENT:
                    ax.plot(t_theory, tsi_replacement_ys, color=common.FIT_COLOR, linewidth=1, label='sDFT fit (t.s.i replace.)' if box_size_index==0 else '')


                if SHOW_PLATEAUS_THEORY:
                    ax.hlines(
                        countoscope_theory.nmsd.plateau_inter_2d(N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)),
                        t[0], t[-1], linestyle='dotted', color=color, label='sDFT plateau' if box_size_index==0 else None)

                if SHOW_PLATEAUS_OBS:
                    plat, plat_std = get_plateau(N2_mean[box_size_index, :], file, L, phi, sigma, t=t_all)
                    ax.hlines(plat, t[0], t[-1], linestyle='dotted', color='grey', label='obs plat' if box_size_index==0 else None)
        
                    if SHOW_PLATEAU_OBS_AREA:
                        ax.fill_between((t[0], t[-1]), [plat-plat_std]*2, [plat+plat_std]*2, facecolor='gray', alpha=0.5)

        assert have_displayed_at_least_one, 'display was false for all L'

        if not RESCALE_X and not FORCE_HIDE_LEGEND:
            ax.legend(fontsize=6 if not PRESENT_SMALL else 6, loc=LEGEND_LOCATION)
        if not LINEAR_Y:
            ax.semilogy()
        ax.semilogx()
        if RESCALE_X == RESCALE_X_L:
            xlabel = '$\Delta t/L$'
        elif RESCALE_X == RESCALE_X_L2:
            xlabel = '$\Delta t/L^2$'
        else:
            xlabel = '$\Delta t$ ($\mathrm{s}$)'
        if RESCALE_Y == False:
            ylabel = r'$\langle \Delta N^2(\Delta t) \rangle$ ($\mathrm{\mu m^2}$)'
        elif RESCALE_Y == RESCALE_Y_N:
            ylabel = r'$\Delta N^2(\Delta t)/ \langle N\rangle$'
        elif RESCALE_Y == RESCALE_Y_VAR_N:
            ylabel = r'$\Delta N^2(\Delta t)/ Var(N)$'
        elif RESCALE_Y == RESCALE_Y_PLATEAU:
            ylabel = r'$\Delta N^2(\Delta t)/ \mathrm{plateau}$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        title = file
        # title = f'Simulated colloids in RCP spheres\n$\phi={phi:.3f}$'
        if not np.isnan(phi):
            title += f', $\phi_\mathrm{{calc}}={phi:.3f}$'
        if not np.isnan(sigma):
            title += f', $\sigma={sigma:.2f}\mathrm{{\mu m}}$'
        if sigma_calced := data.get('particle_diameter_calced'):
            title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
            # print('sigma calced hidden from legend')
        if not PRESENT_SMALL:
            ax.set_title(title)

        filename = f'nmsd_'
        # if SHOW_JUST_ONE_BOX:
        #     filename += f'one_'
        # if RESCALE_X_L2 or RESCALE_Y:
        #     filename += 'rescaled_'
        # if SHOW_T_SLOPE:
        #     filename += f't_'
        # if SHOW_THEORY_FIT:
        #     filename += f'theory_'
        # if SHOW_SHORT_TIME_FIT:
        #     filename += f'shorttime_'
        # if SHOW_TIMESCALEINT_REPLACEMENT:
        #     filename += f'timescaleintreplace_'
        filename += f'{file}'
        if LINEAR_Y:
            filename += '_liny'

        # common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/{filename}.pdf', hide_metadata=True)
        common.save_fig(fig, f'box_counting/figures_png/{filename}.png', dpi=200)

        common.save_data(f'visualisation/data/Ds_from_boxcounting_{file}',
                Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, Ls=Ls_for_saving,
                particle_diameter=sigma,
                pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))
        common.save_data(f'visualisation/data/Ds_from_boxcounting_shorttime_{file}',
                Ds=Ds_shorttime_for_saving, D_uncs=D_uncs_shorttime_for_saving, Ls=Ls_shorttime_for_saving,
                particle_diameter=sigma,
                pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))
        common.save_data(f'visualisation/data/Ds_from_boxcounting_first_quad_{file}',
                Ds=Ds_first_quad_for_saving, D_uncs=D_uncs_first_quad_for_saving, Ls=Ls_first_quad_for_saving,
                particle_diameter=sigma,
                pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))
        common.save_data(f'visualisation/data/Ds_from_boxcounting_collective_{file}',
                Ds=Ds_for_saving_collective, D_uncs=D_uncs_for_saving_collective, Ls=Ls_for_saving_collective,
                particle_diameter=sigma,
                pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))
        
common.save_fig(gradfig, 'box_counting/figures_png/grad.png')