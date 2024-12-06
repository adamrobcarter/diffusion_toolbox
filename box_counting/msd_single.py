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

# if SHOW_JUST_ONE_BOX:
    # LABELS_ON_PLOT = False

# RESCALE_X = RESCALE_X_L2
# RESCALE_X = RESCALE_X_L
RESCALE_X = None

# RESCALE_Y = RESCALE_Y_VAR_N
# RESCALE_Y = RESCALE_Y_PLATEAU
# RESCALE_Y = RESCALE_Y_N
# RESCALE_Y = RESCALE_Y_L2
RESCALE_Y = None

SHORTTIME_FIT_ERROR_THRESH = 0.05 # D_unc/D must be smaller than this for the point to save

have_displayed_at_least_one = False

def get_plateau_range(file):
    start_index = -2500
    end_index   = -1000

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
    elif 'sim_nohydro' in file:
        start_index = -20000
        end_index   = -10000
    else: # used to be -300, -100
        start_index = -600
        end_index   = -400

    return start_index, end_index

gradfig, (gradax1, gradax2, gradax3, gradax4) = plt.subplots(4, 1, figsize=(5, 4*4))
gradax1.semilogx()
gradax2.loglog()
gradax3.loglog()
gradax4.loglog()

# def get_plateau(method, nmsd, file, L, phi, sigma, t, display=False):

#     if display:
#         assert nmsd.size == t.size, f'nmsd {nmsd.shape}, t {t.shape}'
#         grad = np.gradient(nmsd, t)
#         gradax1.scatter(t, grad, s=1)
#         gradax2.scatter(t, np.abs(grad), s=1)

#         t_indexes = np.arange(0, t.size, 1)

#         width = 0.05
#         filter_size = [int(s*width) for s in t_indexes]
#         # print(filter_size[:100])
#         # print(filter_size[::500])
#         # nmsd2 = [np.abs(grad)[i-filter_size[i]:min(i+filter_size[i]+1, t.size-1)].mean() for i in t_indexes]
#         # gradax3.scatter(t, nmsd2, s=1)

#         points = common.exponential_integers(1, t.size-1, 100)
#         # gradax4.scatter(t[points], np.abs(grad)[points], s=1)
#         nmsd3 = [nmsd[points[i]:points[i+1]].mean() for i in range(points.size-1)]
#         # if you change grad to MSD it looks clean AF...
#         gradax4.scatter(t[points][:-1], nmsd3, s=1)

#         # gradax3.scatter(t, scipy.signal.savgol_filter(np.abs(grad)[1000:], 1000, 3))

PLATEAU_SOURCES = ['var', 'varmod', 'nmsdfit', 'nmsdfitinter', 'sDFT', 'N_mean',
                   'density', 'var_time', 'cutoff']
PLATEAU_SOURCE_NAMES = {
    'var': 'variance',
    'varmod': 'variance (single box, mean over boxes)',
    'nmsdfit': 'fit',
    'sDFT': 'theory',
}

def get_plateau(method, nmsd, file, L, phi, sigma, t, var, varmod, D0, N_mean, density, var_time, cutoff_L, cutoff_plat, display=False):
    if method == 'obs':
        return get_plateau_obs(file, nmsd)
    elif method == 'var':
        return 2*var, 0
    elif method == 'varmod':
        return 2*varmod, 0
    elif method == 'var_time':
        return 2*varmod, 0
    elif method == 'nmsdfit':
        return get_plateau_nmsd_fit(nmsd, t, var, L)
    elif method == 'nmsdfitinter':
        return get_plateau_fit_nmsd_inter(file, t, nmsd, phi, sigma, L, D0)
    elif method == 'sDFT':
        return countoscope_theory.nmsd.plateau_inter_2d(N_mean, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)), 0
    elif method == 'N_mean':
        return 2*N_mean, 0
    elif method == 'density':
        return 2*density * L**2, 0
    elif method == 'cutoff':
        if L > cutoff_L:
            return cutoff_plat * L**2 / cutoff_L**2, 0
        else:
            return var*2, 0
    elif method == 'target_nofit':
        pass
    elif method == 'target_fixexponent':
        theory = countoscope_theory.nmsd.plateau_inter_2d(N_mean, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))

        def D_from_plat(plat):
            T_integrand = box_counting.D_of_L.T_integrand_func(nmsd, plat)
            assert np.isfinite(T_integrand).all()
            M1_index, M2_index, _ = box_counting.D_of_L.get_M1_M2(file, L, t, T_integrand)
            _, _, D_of_L, _, _ = box_counting.D_of_L.timescaleintegral_fixexponent(t, L, T_integrand, T_integrand, T_integrand, M2_index, M1_index)
            return D_of_L
        
        # to_minimise = lambda plat: (theory - D_from_plat(plat))**2
        # res = scipy.optimize.minimize(to_minimise, x0=[var*2])
        # return res.x, 0

        plats = np.logspace(np.log10(var/100), np.log10(N_mean), num=2000)
        Ds = np.full_like(plats, np.nan)
        for i in range(len(plats)):
            Ds[i] = D_from_plat(plats[i])

        target_index = np.argmin((Ds-theory)**2)
        # print((Ds-theory)**2)
        # assert target_index != 0
        if target_index == 0:
            warnings.warn('minimisation failed')
            return np.nan, 0

        return plats[target_index], 0


    # elif method == 'S(k=0)<N>':
    #     return
    else:
        assert False



def get_plateau_obs(file, nmsd):
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

    return obs_plat_mean, obs_plat_std

def get_plateau_fit_nmsd_inter(file, t, nmsd, phi, sigma, L, D0):
    # import traceback
    # traceback.print_stack()
        
    obs_plat_mean, _ = get_plateau_obs(file, nmsd) # used for p0 for the optimisation
    
    nmsd_th = countoscope_theory.nmsd.inter_2d(t, D0, 1, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
    nmsd_th_plat = countoscope_theory.nmsd.plateau_inter_2d(1, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
    nmsd_th = nmsd_th / nmsd_th_plat
    
    # points = common.exponential_integers(1, t.size-1, t.size//10)
    points = common.exponential_indices(t, num=t.size//10)
    func = lambda plat : np.sum((plat*nmsd_th[points] - nmsd[points])**2)
    # popt, pcov = scipy.optimize.curve_fit(func, [int(i) for i in np.arange(0, t.size)], nmsd)
    res = scipy.optimize.minimize(func, x0=[obs_plat_mean])
    if not res.success:
        print(res.message)
        if res.message != 'Desired error not necessarily achieved due to precision loss.':
            raise Exception(res.message)
        # assert res.x
    plat_optimised = res.x[0]

    return plat_optimised, 0


def get_plateau_nmsd_fit(nmsd, t, var, L):
    N2_theory2 = lambda t, D, plateau : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
    log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
    
    fitting_points23 = common.exponential_indices(t[1:], num=500)
    # p0 = (0.05, N_mean[box_size_index])
    p02 = [0.05, 2*var]
    popt2, pcov2 = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points23], np.log(nmsd[fitting_points23]), p0=p02, maxfev=2000)
    # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
    # D_from_fit2 = popt2[0]
    # D_from_fit_unc2 = np.sqrt(pcov2[0][0])

    # fit_ys = N2_theory2(t_theory, *popt2)

    return popt2[1], np.sqrt(pcov2[1, 1])

def go(file, ax, separation_in_label=False,
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

    LOWTIME_FIT_END = 20
    
    if rescale_x and rescale_y:
        labels_on_plot = False
        print('I disabled labels on plot as you are rescaling both axes')

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

    num_timesteps = N2_mean.shape[1]
    t_all = data.get('t', np.arange(0, num_timesteps) * time_step)

    D_MSD, _, _ = visualisation.Ds_overlapped.get_D0(file)

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

    plotted_handles = []

    for box_size_index in tqdm.tqdm(iter, desc='box sizes'):
        L   = box_sizes[box_size_index]
        sep = sep_sizes[box_size_index]
        
        L_over_sigma_str = f' = {L/sigma:.2f}σ' if sigma else ''
        print(f'L = {L:.2g}{L_over_sigma_str}')


        delta_N_sq     = N2_mean[box_size_index, :]
        delta_N_sq_err = N2_std [box_size_index, :]
        t = np.copy(t_all)
        t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()*1), 100)


        anomalous = delta_N_sq < 1e-14
        # anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        # print(anomalous)
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
        nmsd_negative = np.sum(delta_N_sq<0)/delta_N_sq.size
        if nmsd_negative: print('nmsd negative', nmsd_negative)
        nmsd_negative = np.sum(delta_N_sq<1e-14)/delta_N_sq.size
        if nmsd_negative: print('nmsd negative', nmsd_negative)
        nmsd_inf = np.sum(np.isinf(delta_N_sq))/delta_N_sq.size
        if nmsd_inf: print('nmsd inf', nmsd_inf)
        
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
        if True: # we calculate the fit even if we don't need to, because we use it for getting the angles for labels_on_plot
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
            
            fitting_points = common.exponential_indices(t, num=500)
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
            # print('  plats', get_plateau(N2_mean[box_size_index, :], file, L, phi, sigma, t=t_all)[0], N2_theory_points[-1])


        
        # # fit to whole thing 2 - replace timescaleint
        # def timescaleint_replacement(nmsd):
        #     # plateau, plateau_unc = get_plateau(plateau_source, nmsd, file, L, phi, sigma, t, var, varmod, display=False)
        #     plateau, plateau_unc = get_plateau(timescaleint_replacement_plateau_source, nmsd, file, L, phi, sigma, t, N_var[box_size_index], N_var_mod[box_size_index], N_mean[box_size_index], D_MSD)
        #     # plateau, plateau_unc = get_plateau_fit_nmsd_inter(file, t, nmsd, phi, sigma, L, D_MSD)
        #     print(f'  plat = {common.format_val_and_unc(plateau, plateau_unc, latex=False)}')
        #     N2_theory2 = lambda t, D : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
        #     log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
            
        #     fitting_points2 = common.exponential_indices(t, num=100)
        #     # p0 = (0.05, N_mean[box_size_index])
        #     p02 = [0.05]
        #     # if file == 'eleanorlong010' and L > 6:
        #     #     p02 = [0.09]
        #     print('dtype', np.log(nmsd).dtype)
        #     print('dtype', np.log(t).dtype)
        #     print('dtype', log_N2_theory2(t, 0.05).dtype)
        #     popt2, pcov2, infodict, mesg, ier = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points2], np.log(nmsd[fitting_points2]), p0=p02, maxfev=2000, full_output=True)
        #     print('mesg', mesg)
        #     print(infodict)
        #     print(ier)
        #     # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
        #     D_from_fit2 = popt2[0]
        #     D_from_fit_unc2 = np.sqrt(pcov2[0][0])

        #     if D_from_fit_unc2/D_from_fit2 < plateau_unc/plateau:
        #         D_from_fit_unc2 = D_from_fit2 * plateau_unc/plateau
        #         # print(f'  increased error to {plateau_unc/plateau:.3f}')

        #     if popt2[0] == p02[0]:
        #         print('  tsi replacement fit failed')
        #         # raise Exception('tsi replacement fit failed')
        #         # D_from_fit2 = np.nan

        #     fit_ys = N2_theory2(t_theory, *popt2)

        #     return D_from_fit2, D_from_fit_unc2, fit_ys
        

        plateau, plateau_unc = get_plateau(
            timescaleint_replacement_plateau_source, 
            nmsd=delta_N_sq,
            file=file,
            L=L,
            phi=phi,
            sigma=sigma,
            t=t,
            var=N_var[box_size_index],
            varmod=N_var_mod[box_size_index],
            N_mean=N_mean[box_size_index],
            D0=D_MSD,
            density=data.get('density', np.nan),
            var_time=None,
            cutoff_L=None,
            cutoff_plat=None
        )
            
        # fitting_points2 = common.exponential_indices(t, num=100)
        # # fitting_points2

        # N2_theory2 = lambda t, D, plateau : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
        # log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
            

        def timescaleint_replacement_fitplateau():
            assert np.isnan(delta_N_sq).sum() == 0, 'nan found in nmsd'
            assert np.sum(delta_N_sq < 0) == 0, 'negatives found in nmsd'
            assert np.sum(delta_N_sq == 0) == 0, 'zeros found in nmsd'

            N2_theory2 = lambda t, D, plateau : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
            log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
            

            t_for_fit = t[1:]
            t_for_fit = t_for_fit[t_for_fit > sigma**2/D_MSD]
            # t_for_fit2 = t_for_fit[t_for_fit > 3*L**2]
            # if t_for_fit2.size > 100:
            #     t_for_fit = t_for_fit2

            fitting_points22 = common.exponential_indices(t_for_fit, num=100)
            # p0 = (0.05, N_mean[box_size_index])
            p02 = [0.05, 2*N_var[box_size_index]]


            
            assert np.isfinite(t[fitting_points22]).all()
            assert np.isfinite(np.log(delta_N_sq[fitting_points22])).all()
            popt2, pcov2 = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points22], np.log(delta_N_sq[fitting_points22]), p0=p02, maxfev=2000)
            # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
            D_from_fit2 = popt2[0]
            D_from_fit_unc2 = np.sqrt(pcov2[0][0])

            fit_ys = N2_theory2(t_theory, *popt2)

            return D_from_fit2, D_from_fit_unc2, fit_ys
        
        def timescaleint_replacement():
            assert np.isnan(delta_N_sq).sum() == 0, 'nan found in nmsd'
            assert np.sum(delta_N_sq < 0) == 0, 'negatives found in nmsd'
            assert np.sum(delta_N_sq == 0) == 0, 'zeros found in nmsd'

            N2_theory2 = lambda t, D : plateau * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
            log_N2_theory2 = lambda t, *args : np.log(N2_theory2(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
            

            t_for_fit = t[1:]
            # t_for_fit = t_for_fit[t_for_fit > sigma**2/D_MSD]
            tmax1 = t_for_fit.max()
            # t_for_fit = t_for_fit[t_for_fit < 1e2]
            tmax2 = t_for_fit.max()
            # assert tmax2 < tmax1
            # print('t fitting 2', t_for_fit[0], t_for_fit[-1])

            # t_for_fit2 = t_for_fit[t_for_fit > 3*L**2]
            # if t_for_fit2.size > 100:
            #     t_for_fit = t_for_fit2

            fitting_points22 = common.exponential_indices(t_for_fit, num=100)
            # p0 = (0.05, N_mean[box_size_index])
            assert t[fitting_points22[-1]] < tmax1
            # print('t fitting 2', t[fitting_points22[0]], t[fitting_points22[-1]])
            p02 = [0.05]


            
            assert np.isfinite(t[fitting_points22]).all()
            assert np.isfinite(np.log(delta_N_sq[fitting_points22])).all()
            popt2, pcov2 = scipy.optimize.curve_fit(log_N2_theory2, t[fitting_points22], np.log(delta_N_sq[fitting_points22]), p0=p02, maxfev=2000)
            # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
            D_from_fit2 = popt2[0]
            D_from_fit_unc2 = np.sqrt(pcov2[0][0])

            fit_ys = N2_theory2(t_theory, *popt2)

            return D_from_fit2, D_from_fit_unc2, fit_ys
        
        # Dc, Dc_unc, tsi_replacement_ys = timescaleint_replacement_fitplateau()
        Dc, Dc_unc, tsi_replacement_ys = timescaleint_replacement()
        
        print(f'  timescaleint replacement D={common.format_val_and_unc(Dc, Dc_unc, latex=False)}')
        # Dc_lower, Dc_unc_lower, _ = tsi_replace_func((delta_N_sq-delta_N_sq_err).clip(min=0.1)) # otherwise we could give negatives
        # Dc_upper, Dc_unc_upper, _ = tsi_replace_func(delta_N_sq+delta_N_sq_err)
        # Dc_unc_final = max(Dc_unc, abs(Dc-Dc_lower), abs(Dc-Dc_upper))
        Dc_unc_final = Dc_unc # errors are way overestimated with the above lines
        # print(f'  final Dc unc {Dc_unc_final/Dc:.3f}')

        Ds_for_saving_collective.append(Dc)
        D_uncs_for_saving_collective.append(Dc_unc_final)
        Ls_for_saving_collective.append(L)
        
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
            N2_theory_points    /= rescale_y_value
            N2_theory_points_Lh /= rescale_y_value
            tsi_replacement_ys  /= rescale_y_value

        if rescale_x:
            rescale_x_value = 1
            if rescale_x == RESCALE_X_L2:
                rescale_x_value = L**2
            elif rescale_x == RESCALE_X_L:
                rescale_x_value = L
            
            t_theory = t_theory / rescale_x_value
            t        = t        / rescale_x_value
        
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
        
        if box_size_indices:
            display = box_size_index in box_size_indices
        elif SHOW_JUST_ONE_BOX:
            display = box_size_index == 20
        elif len(box_sizes) <= max_boxes_on_plot:
            display = True
        else:
            display = box_size_index % (len(box_sizes) // max_boxes_on_plot) == 0

        if display:
            print('box_size_index', box_size_index)
            # get_plateau(nmsd=N2_mean[box_size_index, :], file=file, L=L, phi=phi, sigma=sigma, display=True, t=t_all)
            
            if LINEAR_Y:
                if display_i != 0:
                    ax = ax.twinx() 
                yrange = delta_N_sq.max()
                ylim_range = yrange * max_boxes_on_plot
                OVERLAP = 0.5
                ylim_start =  - OVERLAP * yrange * display_i
                ylim_end = yrange * (1 + (max_boxes_on_plot - display_i)*OVERLAP)
                ax.set_ylim(ylim_start, ylim_end)
                ax.get_yaxis().set_visible(False)
            
            display_i += 1

            have_displayed_at_least_one = True

            if SHOW_MEAN:
                ax.hlines(2*N_mean[box_size_index]/rescale_y_value, t.min(), t.max(), color=color, linewidth=1, linestyle='dashdot', label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
            if show_variance:
                ax.hlines(2*N_var [box_size_index]/rescale_y_value, t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$2\mathrm{Var}(N)$' if box_size_index==0 else None)
                ax.hlines(2*N_var_mod[box_size_index]/rescale_y_value, t.min(), t.max(), linestyles='dotted', color='grey', linewidth=1, label=r'$2\mathrm{Var*}(N)$' if box_size_index==0 else None)


            if not rescale_y and SHOW_THEORY_FIT:
                ax.plot(t_theory[1:], N2_theory_points[1:], color='gray', linewidth=1, label=type_of_fit if box_size_index==0 else None)
                ax.plot(t_theory[1:], N2_theory_points_Lh[1:], color='white', linewidth=1, label=type_of_fit if box_size_index==0 else None)
            
            if show_rescaled_theory and box_size_index==0:
                t_theory_rescaled = np.logspace(-4, 4, num=200)
                theory = N2_theory(t_theory_rescaled, D_MSD)
                ax.plot(t_theory_rescaled[1:]/rescale_x_value, theory[1:]/rescale_y_value, color='black', linewidth=1, label='theory')
                
            
            # label += fr', $D_\mathrm{{fit}}={popt[0]:.3g}\pm {np.sqrt(pcov[0][0]):.3g} \mathrm{{\mu m^2/s}}$'

            if (rescale_x or rescale_y): # remove and False in future please
                markersize = 2
            else:
                if PRESENT_SMALL:
                    markersize = 5
                else:
                    markersize = 3

            # actual data
            if labels_on_plot:
                label = label='observations' if box_size_index==0 else ''
            else:
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
            print('  dN2s', delta_N_sq.size, common.nanfrac(delta_N_sq))

            if DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE and delta_N_sq.size > 1000:
                # print('tt', t)
                points_to_plot = common.exponential_indices(t, 500)
            else:
                points_to_plot = np.index_exp[1:] # this is basically a functional way of writing points_to_plot = [1:]
            
            exp_plot, = ax.plot(t[points_to_plot], delta_N_sq[points_to_plot], label=label, linestyle=linestyle, marker='o', markersize=markersize, zorder=-1, color=color)
            # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:]/np.sqrt(num_of_boxes[box_size_index]), label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            plotted_handles.append(exp_plot)

            if labels_on_plot:
                t_index_for_text = int(t_theory.size // 1.6)
                angle = np.arctan(np.gradient(N2_theory_points, t_theory)[t_index_for_text]) * 180/np.pi
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

            if show_timescaleint_replacement:
                ax.plot(t_theory, tsi_replacement_ys, color=common.FIT_COLOR, linewidth=1, label='sDFT fit' if box_size_index==0 else '')


            if SHOW_PLATEAUS_THEORY:
                ax.hlines(
                    countoscope_theory.nmsd.plateau_inter_2d(N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)),
                    t[0], t[-1], linestyle='dotted', color=color, label='sDFT plateau' if box_size_index==0 else None)

            if SHOW_PLATEAUS_OBS:
                plat, plat_std = get_plateau(nmsd=N2_mean[box_size_index, :], file=file, L=L, phi=phi, sigma=sigma, D0=D_MSD, t=t_all)
                ax.hlines(plat, t[0], t[-1], linestyle='dotted', color='grey', label='obs plat' if box_size_index==0 else None)
    
                if SHOW_PLATEAU_OBS_AREA:
                    ax.fill_between((t[0], t[-1]), [plat-plat_std]*2, [plat+plat_std]*2, facecolor='gray', alpha=0.5)

    if show_nointer_theory_limits:
        assert rescale_x == RESCALE_X_L2
        t_theory_limits_over_L2 = np.logspace(-2, 3)
        L_low = 0.01
        print('pack frac', data['pack_frac'])
        Dc = D_MSD * (1 + data['pack_frac']) / (1 - data['pack_frac'])**3
        print('Dc/D_MSD', Dc/D_MSD)
        L_high = 1000
        handles = []
        for i, L, D, linest in ((0, L_low, D_MSD, 'dotted'), (1, L_high, Dc, 'dashed')):
        # for L, D in ((L_high, Dc),):
            t_theory_limits = t_theory_limits_over_L2 * L**2
            # N_mean_low = data['density'] * L_low
            N_mean = 0.01643333546467505 * L**2
            plateau_over_2 = 1/2 * countoscope_theory.nmsd.plateau_inter_2d(N_mean, L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
            print('we', D, N_mean, L_low)
            nmsd_limit = countoscope_theory.nmsd.nointer_2d(t_theory_limits, D, plateau_over_2, L)
            line, = ax.plot(t_theory_limits_over_L2[1:], nmsd_limit[1:]/L**2, color='black', linewidth=1, linestyle=linest)#, label=nointer_theory_limit_labels[i])                
            print('t', t_theory_limits_over_L2[1:])
            print('N', nmsd_limit[1:])
            handles.append(line)
        # legend2 = ax.legend(handles, labels=['a', 'b'])
        print(handles)
        if show_second_legend:
            legend2 = ax.legend(handles=handles, labels=nointer_theory_limit_labels, loc='lower right', bbox_to_anchor=(0, 0, 0.75, 1), fontsize=legend_fontsize)
            ax.add_artist(legend2)

    assert have_displayed_at_least_one, 'display was false for all L'

    if rescale_x and not FORCE_HIDE_LEGEND:
        legend = ax.legend(handles=plotted_handles, fontsize=legend_fontsize, loc=legend_location)
        common.set_legend_handle_size(legend)
        pass
    if not LINEAR_Y:
        ax.semilogy()
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
    # filename += f'{file}'
    # if LINEAR_Y:
    #     filename += '_liny'

    # common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/{filename}.pdf', hide_metadata=True)
    # common.save_fig(fig, f'box_counting/figures_png/{filename}.png', dpi=200)
    # if export_destination:
    #     common.save_fig(fig, export_destination, hide_metadata=True)

    common.save_data(f'visualisation/data/Ds_from_boxcounting_{file}',
            Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, Ls=Ls_for_saving,
            particle_diameter=sigma,
            pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'), max_time_hours=data.get('max_time_hours'))
    common.save_data(f'visualisation/data/Ds_from_boxcounting_shorttime_{file}',
            Ds=Ds_shorttime_for_saving, D_uncs=D_uncs_shorttime_for_saving, Ls=Ls_shorttime_for_saving,
            particle_diameter=sigma,
            pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'), max_time_hours=data.get('max_time_hours'))
    common.save_data(f'visualisation/data/Ds_from_boxcounting_first_quad_{file}',
            Ds=Ds_first_quad_for_saving, D_uncs=D_uncs_first_quad_for_saving, Ls=Ls_first_quad_for_saving,
            particle_diameter=sigma,
            pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'), max_time_hours=data.get('max_time_hours'))
    common.save_data(f'visualisation/data/Ds_from_boxcounting_collective_{timescaleint_replacement_plateau_source}_{file}',
            Ds=Ds_for_saving_collective, D_uncs=D_uncs_for_saving_collective, Ls=Ls_for_saving_collective,
            particle_diameter=sigma,
            pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'), max_time_hours=data.get('max_time_hours'))
        
    print('largest box', box_sizes[-1])
# common.save_fig(gradfig, 'box_counting/figures_png/grad.png')

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
           show_title=not PRESENT_SMALL,
           labels_on_plot=LABELS_ON_PLOT,
           rescale_x=RESCALE_X, rescale_y=RESCALE_Y,
           max_boxes_on_plot=MAX_BOXES_ON_PLOT,
        )

        ax.set_ylim(1e-3, 5e2)
        if RESCALE_X == RESCALE_X_L2 and RESCALE_Y == RESCALE_Y_L2:
            ax.set_ylim(3e-3, 0.7e-1)
            ax.set_xlim(1e-2, 1e4)
        
        common.save_fig(fig, f'box_counting/figures_png/nmsd_{file}.png', dpi=200)