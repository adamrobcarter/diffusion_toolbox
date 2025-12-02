import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate, scipy.special, scipy.optimize, scipy.signal
import countoscope_theory.nmsd, countoscope_theory.structure_factor
import tqdm
import math
import matplotlib.cm
import box_counting.msd_single
import visualisation.Ds_overlapped

SHOW_THEORY = True
SHOW_TIMESCALEINTEGRAL_FIT_SHORT = False
SHOW_TIMESCALEINTEGRAL_FIT_LONG = True
SHOW_THRESH_LINE = True
SHOW_NOFIT_CUTOFF = True
DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE = True
SHOW_LEGEND = True

LATE_CN_ALPHA = 1.0
LABELS_ON_PLOT = False
LABELS_ON_PLOT_Y_SHIFT = 1.2
LABELS_ON_PLOT_X_SHIFT = 1.3

MAX_NUM_BOXES = 10

UNC_INCLUDES_NMSD_UNC = False
UNC_INCLUDES_VAR_UNC  = True

NOFIT_CROP_THRESHOLD = 1e-4

RESCALE_X_L2 = 1
# RESCALE_X = RESCALE_X_L2
RESCALE_X = None

# PLATEAU_SOURCE = 'cutoff'
# PLATEAU_SOURCE = 'histogram'
PLATEAU_SOURCE = 'var'
# PLATEAU_SOURCE = 'varmod'
# PLATEAU_SOURCE = 'nmsdfit'
# PLATEAU_SOURCE = 'sDFT'
# PLATEAU_SOURCE = 'target_fixexponent'
# PLATEAU_SOURCE = 'nmsdfitinter'
# PLATEAU_SOURCE = 'var_losecorr'

"""
Countoscope appendix:

Given the above estimates for h N2(t)i and hN2i hNi2, the timescale integral in equation (8) 
is calculated using trapezoidal summation. When the observation box size is very small, the 
integrand in Eq. (8) decays very quickly, and much of its support may fall before the first 
frame. When the observation box is very large, the decay is slow and we may not have enough 
data to calculate the integral accurately. Both of these issues can be addressed by fitting 
our available data at short and long times to theory-informed functional forms and calculating 
the missing contributions to the timescale integral. This procedure is described in detail in 
Supp. Mat. Sec. 1.5 which also includes a discussion of how the errorbars were calculated for 
Fig. 5.
"""
f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2

import warnings
# warnings.filterwarnings('ignore')

def do_early_integral(f, t, T_integrand, T_integrand_min, T_integrand_max):
    interval_end_index = 4
    early_fit_func = lambda t, a: f(t/a)**4
    early_fit_xdata = [0, t[1]]
    early_fit_ydata = [1, T_integrand[1]]
    assert np.isfinite(early_fit_xdata[1])
    assert np.isfinite(early_fit_ydata[1])

    early_popt, early_pcov = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata)
    early_integral = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt[0])[0]
        
    early_fit_ydata_min = [1, T_integrand_min[1]]
    assert np.isfinite(early_fit_ydata_min[1])
    early_popt_min, early_pcov_min = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata_min)
    early_integral_oneway = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt_min[0])[0]
        
    early_fit_ydata_max = [1, T_integrand_max[1]]
    assert np.isfinite(early_fit_ydata_max[1])
    early_popt_max, early_pcov_max = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata_max)
    early_integral_otherway = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt_max[0])[0]
    
    # if UNC_INCLUDES_NMSD_UNC:
    early_integral_min = min(early_integral, early_integral_oneway, early_integral_otherway)
    early_integral_max = max(early_integral, early_integral_oneway, early_integral_otherway)
    # else:
    #     early_integral_min = early_integral
    #     early_integral_max = early_integral

    assert early_integral_min <= early_integral
    assert early_integral <= early_integral_max
    
    return early_fit_func,early_popt,early_integral,early_integral_min,early_integral_max

class PositiveSlopeInTimescaleIntegralFitException(Exception):
    pass

def do_late_integral(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index, fix_slope=False):
    if fix_slope:
        b = 4
        late_fit_func = lambda t, a, B: (a - b*t) - np.log(1 + np.exp(a - b*t))
    else:
        late_fit_func = lambda t, a, b: (a - b*t) - np.log(1 + np.exp(a - b*t))
    late_fit_xdata =     np.log(t          [M1_index:M2_index])
    late_fit_ydata = 2 * np.log(T_integrand[M1_index:M2_index])
    assert late_fit_xdata.size > 3
    late_popt, late_pcov = scipy.optimize.curve_fit(late_fit_func, late_fit_xdata, late_fit_ydata)
    if fix_slope:
        a = late_popt[0]
    else:
        a, b = late_popt
    # print(f'a = {a}, b = {b}')
    if b < 0:
        print('  positive slope in timescale integral fit')
        raise PositiveSlopeInTimescaleIntegralFitException()
    a_unc = np.sqrt(late_pcov)[0, 0]
    b_unc = np.sqrt(late_pcov)[1, 1]
    late_fit_func_real = lambda t, a, b: np.exp(0.5 * late_fit_func(np.log(t), a, b))
    # what is real? is the other one log space or something?
    t_end = t.max()*1e3
    if fix_slope:
        t_fit_late = np.logspace(np.log10(t[M2_index]), np.log10(t_end))
    else:
        t_fit_late = np.logspace(np.log10(t[M1_index]), np.log10(t_end))
    if late_fit_func_real(t_end, a, b) < 1e-5:
        warnings.warn(f'late fit has not decayed to zero')# = late_fit_func_real({t_end}, {a}, {b}) = {late_fit_func_real(t_end, a, b)}')
    
    assert np.all(late_fit_func_real(t_fit_late, *late_popt) > 0), f'some of late_fit_func_real were negative for L={L}'
    assert np.all(late_fit_func_real([t[M2_index], t_end], *late_popt) > 0), f'some of late_fit_func_real were negative for L={L}'
    assert np.all(late_fit_func_real(np.logspace(np.log10(t[M2_index]), np.log10(t_end), 1000), *late_popt) > 0), f'some of late_fit_func_real were negative for L={L}'
    assert np.all(late_fit_func_real(np.linspace(t[M2_index], t_end, 1000), *late_popt) > 0), f'some of late_fit_func_real were negative for L={L}'
    assert t[M2_index] < t_end
    p     = scipy.integrate.quad(late_fit_func_real, t[M2_index], t_end, args=(a, b), full_output=True)
    late_integral = p[0]
    # assert np.all(late_fit_func_real(p[2]['alist'], *late_popt) > 0)
    # assert np.all(late_fit_func_real(p[2]['blist'], *late_popt) > 0)
    late_integral_00 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t_end, args=(a+a_unc, b+b_unc))[0]
    late_integral_01 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t_end, args=(a+a_unc, b-b_unc))[0]
    late_integral_10 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t_end, args=(a-a_unc, b+b_unc))[0]
    late_integral_11 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t_end, args=(a-a_unc, b-b_unc))[0]

    int_min = min(late_integral, late_integral_00, late_integral_01, late_integral_10, late_integral_11)
    int_max = max(late_integral, late_integral_00, late_integral_01, late_integral_10, late_integral_11)
    if late_integral != 0:
        # print(f'  late integral uncs {int_min/late_integral:.3f} {int_max/late_integral:.3f}')
        pass
    else:
        warnings.warn('late integral = 0')
    t_to_inf = np.logspace(np.log10(t[M1_index]), 15, 1000)
    lf = late_fit_func_real(t_to_inf, *late_popt)
    assert np.all(lf >= 0)

    if late_integral < 0:
        late_integral = 0
        warnings.warn('TURNED OFF NEGATIVE INTEGRAL')
    if int_min < 0:
        int_min = 0
        warnings.warn('TURNED OFF NEGATIVE INTEGRAL')
    if int_max < 0:
        int_max = 0
        warnings.warn('TURNED OFF NEGATIVE INTEGRAL')
    assert late_integral >= 0, f'late_integral = {late_integral}, late_fit_func_real(t[M2_index]) = {late_fit_func_real(t[M2_index], *late_popt)}, late_fit_func_real(t_end) = {late_fit_func_real(t_end, *late_popt)}'

    return late_integral, int_min, int_max, late_popt, late_fit_func_real, t_fit_late

def do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, M2_index):
    assert t.size == T_integrand.size == T_integrand_min.size == T_integrand_max.size

    data_integral = scipy.integrate.trapezoid(T_integrand    [1:M2_index], t[1:M2_index])
    min_and_max = [
            scipy.integrate.trapezoid(T_integrand    [1:M2_index], t[1:M2_index]),
            scipy.integrate.trapezoid(T_integrand_max[1:M2_index], t[1:M2_index]),
            scipy.integrate.trapezoid(T_integrand_min[1:M2_index], t[1:M2_index])
        ]
    # if UNC_INCLUDES_NMSD_UNC:
    data_integral_min = min(min_and_max)
    data_integral_max = max(min_and_max)
    # else:
    #     data_integral_min = data_integral
    #     data_integral_max = data_integral
    assert data_integral_min <= data_integral
    assert data_integral_max >= data_integral
    if data_integral_min == data_integral:
        warnings.warn('slightly strange uncertainty behavoir')
    return data_integral,data_integral_min,data_integral_max

def get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral, data_integral_min, data_integral_max, late_integral, late_integral_min, late_integral_max):
    assert early_integral_min <= early_integral <= early_integral_max
    assert data_integral_min <= data_integral <= data_integral_max
    assert late_integral_min <= late_integral <= late_integral_max

    assert early_integral     >= 0
    assert early_integral_min >= 0
    assert early_integral_max >= 0
    assert data_integral      >= 0
    assert data_integral_min  >= 0
    assert data_integral_max  >= 0
    assert late_integral      >= 0
    assert late_integral_min  >= 0
    assert late_integral_max  >= 0

    total_integral     = early_integral     + data_integral     + late_integral
    total_integral_min = early_integral_min + data_integral_min + late_integral_min
    total_integral_max = early_integral_max + data_integral_max + late_integral_max
    assert total_integral_min <= total_integral
    assert total_integral_max >= total_integral

    T_of_L     = 2 * (total_integral) # countoscope eq. 8
    T_of_L_min = 2 * total_integral_min
    T_of_L_max = 2 * total_integral_max
    assert T_of_L_min <= T_of_L
    assert T_of_L_max >= T_of_L

    alpha_T = 0.561244
    D_of_L     = alpha_T * L**2 / (4 * T_of_L) # countoscope eq. 9
    D_of_L_min = alpha_T * L**2 / (4 * T_of_L_max)
    D_of_L_max = alpha_T * L**2 / (4 * T_of_L_min)

    assert D_of_L_min <= D_of_L, f'D_of_L_min = {D_of_L_min} !<= D_of_L = {D_of_L}'
    assert D_of_L_max >= D_of_L, f'D_of_L_max = {D_of_L_max} !>= D_of_L = {D_of_L}'
    return D_of_L, D_of_L_min, D_of_L_max


def C_N_simplefit(t, C_N_over_VarN, C_N_over_VarN_unc, L):
    
    func = lambda t, D : countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2

    fitting_points = common.exponential_indices(t, num=100)
    p0 = [0.05]
    popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], C_N_over_VarN[fitting_points], sigma=C_N_over_VarN_unc[fitting_points], p0=p0, maxfev=2000)
    
    D_from_fit = popt[0]
    D_from_fit_unc = np.sqrt(pcov[0][0])
    return D_from_fit, D_from_fit_unc, func

def T_integrand_func(nmsd, nmsd_unc, plat, plat_unc):
    T_integrand = (1 - nmsd / plat )**2 # countoscope paper eq. 8
    # T_integrand_unc_sq_nmsd = (2*(1 - nmsd/plat) * (nmsd/plat) * plat_unc)**2
    # T_integrand_unc_sq_var  = (2*(1-nmsd/plat) * (-1/plat) * nmsd_unc)**2       # these I think are wrong but they do look good
    T_integrand_unc_sq_nmsd = (2*(1 - nmsd/plat) * (-1/plat) * nmsd_unc)**2
    T_integrand_unc_sq_var  = (2*(1 - nmsd/plat) * (nmsd/plat**2) * plat_unc)**2

    T_integrand_unc_sq = np.zeros_like(T_integrand)
    if UNC_INCLUDES_NMSD_UNC:
        T_integrand_unc_sq += T_integrand_unc_sq_nmsd
    if UNC_INCLUDES_VAR_UNC:
        T_integrand_unc_sq += T_integrand_unc_sq_var

    return T_integrand, np.sqrt(T_integrand_unc_sq)

def get_M1_M2(file, L, t, T_integrand, plateau_source):
        # data integral setup
    MIN_M2_INDEX = 10

    M2_index = None
    M1_index = None

    thresh_line = None

    if file.startswith('eleanorlong') or file.startswith('eleanor0.01') or file.startswith('brennan'):
        M2_index = min(int(300 * L), t.shape[0]-1)
        if file == 'eleanorlong066':
            M2_index = min(int(2400 * np.sqrt(L)), t.shape[0]-1)
        if file == 'eleanorlong001':
            if plateau_source == 'var':
                M2_index = min(int(5000 * np.sqrt(L)), t.shape[0]-1)
            elif plateau_source == 'nmsdfitinter':
                M2_index = min(int(4000 * np.sqrt(L)), t.shape[0]-1)
            elif plateau_source == 'nmsdfit':
                M2_index = min(int(4000 * np.sqrt(L)), t.shape[0]-1)
            else:
                M2_index = min(int(2600 * np.sqrt(L)), t.shape[0]-1)
            print('M2_index', M2_index)
        M2_index = max(M2_index, MIN_M2_INDEX)
        M1_index = int(M2_index / 1.2) # t is linear so we can just halve the index to halve the time

        if 'eleanorlong010' in file:
            c = 4e-6
            m = 1
            thresh_line = c * t**m
            M2_index = np.argmax(T_integrand < thresh_line)
            # if M2_index < 100:
            #    M2_index = min(int(60 * L), t.shape[0]-1)
            M2_index = max(M2_index, MIN_M2_INDEX)
            assert M2_index != 0
            M1_t = t[M2_index] / 2
            # print('t', t, 'M1_t', M1_t)
            M1_index = np.argmax(t > M1_t)
        
    elif file == 'alice0.02':
        M2_index = min(int(60 * L), t.shape[0]-1)
        M2_index = max(M2_index, MIN_M2_INDEX)
        M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        # D0 = 0.0416
    elif file == 'alice0.02_overlapped3':
        M2_index = min(int(100 * L), t.shape[0]-1)
        M2_index = max(M2_index, MIN_M2_INDEX)
        M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        # D0 = 0.0416
    elif file == 'alice0.34':
        M2_index = min(int(40 * L), t.shape[0]-1)
        M2_index = max(M2_index, MIN_M2_INDEX)
        M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        # D0 = 0.0310
    elif file == 'alice0.66':
        M2_index = min(int(200 * L), t.shape[0]-1)
        M2_index = max(M2_index, MIN_M2_INDEX)
        M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        # D0 = 0.0175
    elif file.startswith('sim_'):
        if 'hydro_010_L640_div8' in file:
            M2_index = min(int(600 * np.sqrt(L)), t.shape[0]-1)
        elif 'hydro_011_L320_longer' in file: # should this be longer?
            # M2_index = min(int(560 * L**0.4), t.shape[0]-1)
            # print(N2_mean[box_size_index, :] < thresh_line)
            c = 4e-6
            m = 1
            thresh_line = c * t**m
            M2_index = np.argmax(T_integrand < thresh_line)
            assert M2_index != 0
            M1_t = t[M2_index] / 2
            M1_index = np.argmax(t > M1_t)
            
        elif 'L320' in file and 'longer' in file and '011' in file:
            c = 1e-6
            m = 1
            thresh_line = c * t**m
            M2_index = np.argmax(T_integrand < thresh_line)
            # assert M2_index != 0
            # print('you removed this assertion which may be bad')
            # if M2_index == 0: # T integrand and thresh_line never crossed
            #     M2_index = t.size-1
            M1_t = t[M2_index] / 2
            M1_index = np.argmax(t > M1_t)
            
        # elif 'hydro_010_L640_longer' in file or 'hydro_011_L640_longer' in file:
        elif 'L640' in file and 'longer' in file and '011' in file:
            # M2_index = min(int(560 * L**0.4), t.shape[0]-1)
            # print(N2_mean[box_size_index, :] < thresh_line)
            c = 5e-7
            m = 1
            thresh_line = c * t**m
            M2_index = np.argmax(T_integrand < thresh_line)
            # assert M2_index != 0
            # print('you removed this assertion which may be bad')
            # if M2_index == 0: # T integrand and thresh_line never crossed
            #     M2_index = t.size-1
            M1_t = t[M2_index] / 2
            M1_index = np.argmax(t > M1_t)
            
        elif 'L1280' in file and 'longer' in file and '011' in file:
            # M2_index = min(int(560 * L**0.4), t.shape[0]-1)
            # print(N2_mean[box_size_index, :] < thresh_line)
            c = 1e-7
            m = 1
            thresh_line = c * t**m
            M2_index = np.argmax(T_integrand < thresh_line)
            assert M2_index != 0
            M1_t = t[M2_index] / 2
            M1_index = np.argmax(t > M1_t)

        elif 'L640' in file and 'longer' in file and '002' in file:
            c = 2e-7
            m = 1
            thresh_line = c * t**m
            M2_index = np.argmax(T_integrand < thresh_line)
            assert M2_index != 0
            M1_t = t[M2_index] / 2
            M1_index = np.argmax(t > M1_t)

        elif 'hydro_002_L320_longer_merged' in file:
            # M2_index = min(int(560 * np.sqrt(L)), t.shape[0]-1)
            M2_index = np.argmax(T_integrand < thresh_line)
            M1_t = t[M2_index] / 2
            M1_index = np.argmax(t > M1_t)
        elif 'hydro_011_L320_longer' in file:
            M2_index = min(int(80 * np.sqrt(L)), t.shape[0]-1)
            M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        else:
            # raise Exception()
            M2_index = min(int(1800 * np.sqrt(L)), t.shape[0]-1)
            M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        M2_index = max(M2_index, MIN_M2_INDEX)

    elif file.startswith('ld_'):
        # c = 1e-7
        # m = 1
        # thresh_line = c * t**m
        # M2_index = np.argmax(T_integrand < thresh_line)
        # assert M2_index != 0
        # M1_t = t[M2_index] / 2
        # M1_index = np.argmax(t > M1_t)
        
        M2_index = min(int(1800 * np.sqrt(L)), t.shape[0]-1)
        M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
        
    else:
        raise Exception('you need to define the parameters for this dataset')

    assert M2_index - M1_index > 3, f'M2_index = {M2_index}, M1_index = {M1_index}, t.shape[0]-1 = {t.shape[0]-1}'

    return M1_index, M2_index, thresh_line

def go(file, plateau_source, ax=None, legend_fontsize=8, title=None, save_data=False,
       plot_color=None, export_destination=None, sep_in_label=False,
       show_nofit_cutoff=False, show_theory=False, labels_on_plot=True, max_num_boxes=10,
       show_short_fits=True, show_long_fits=True, late_C_N_alpha=LATE_CN_ALPHA, labels_on_plot_font_color=None,
       show_legend=False, plot_C_N_squared=True, box_size_indices=None,
       show_slope=False, show_T_of_L=False, rescale_x=None, colormap=None,
       disable_ylabel=False, plateau_adjust=1.0, plateau_source_suffix='', plateau_offset=0.0,
       show_thresh_line=False, hide_data_after_nofit_cutoff=True):

    # ax.set_title(f'{file}: plateau: {plateau_source}')

    # D_fig = plt.figure()
    # T_fig = plt.figure()
    # D_ax = D_fig.gca()
    # T_ax = T_fig.gca()

    # integrand_fig = plt.figure()
    # ax = integrand_fig.gca()
    # integrand_fig, (ax_old, mid_ax, extra_ax) = plt.subplots(3, 1, figsize=(6, 15))
    

    # integrand_rescaled_fig = plt.figure()
    # integrand_rescaled_ax = integrand_rescaled_fig.gca()

    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean   = data['N2_mean']
    N2_std    = data['N2_std']
    sigma     = data['particle_diameter']
    phi       = data['pack_frac']
    time_step = data['time_step']

    box_sizes = data['box_sizes']
    N_mean    = data['N_mean']
    N_var     = data['N_var']
    sep_sizes = data['sep_sizes']

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t = data.get('t', np.arange(0, num_timesteps) * time_step)

    D_of_Ls     = np.full((num_boxes), np.nan)
    D_of_L_uncs = np.full((2, num_boxes), np.nan)
    T_of_Ls     = np.full((num_boxes), np.nan)

    D_of_Ls_fixexponent     = np.full((num_boxes), np.nan)
    D_of_L_uncs_fixexponent = np.full((2, num_boxes), np.nan)
    
    D_of_Ls_nofit     = np.full((num_boxes), np.nan)
    D_of_L_uncs_nofit = np.full((2, num_boxes), np.nan)
    
    D_of_Ls_nofit2     = np.full((num_boxes), np.nan)
    D_of_L_uncs_nofit2 = np.full((2, num_boxes), np.nan)
    
    D_of_Ls_simplefit     = np.full((num_boxes), np.nan)
    D_of_L_uncs_simplefit = np.full((2, num_boxes), np.nan)
    
    D_of_Ls_nofit_noshortfit     = np.full((num_boxes), np.nan)
    D_of_L_uncs_nofit_noshortfit = np.full((2, num_boxes), np.nan)
    
    D0, _, _ = visualisation.Ds_overlapped.get_D0(file)

    for box_size_index in tqdm.trange(num_boxes, desc='box sizes'):
        L = box_sizes[box_size_index]
        
        if sigma:
            L_label = f'L={L/sigma:.2g}σ'
        else:
            L_label = f'L={L:.2g}'
        print(f'{L_label} N_var={N_var[box_size_index]:.4f}, N_var/L^2={N_var[box_size_index]/L**2:.4f}')

        # fits etc
        cutoff_index = np.argmax(box_sizes > 30)
        plateau, plateau_std = box_counting.msd_single.get_plateau(
            plateau_source,
            file=file,
            data=data,
            box_size_index=box_size_index,
            # nmsd=N2_mean[box_size_index, :],
            # L=L,
            phi=phi,
            sigma=sigma,
            t=t,
            # var=N_var[box_size_index],
            # var_std=N_var_std[box_size_index],
            # var_std=N_var[box_size_index] - N_var_sem_lb[box_size_index],
            # var_std=N_var_sem[box_size_index],
            # varmod=N_var_mod[box_size_index],
            D0=D0,
            # N_mean=N_mean[box_size_index],
            # density=data.get('density', np.nan),
            # var_time=data.get('var_time'),
            cutoff_L=box_sizes[cutoff_index],
            cutoff_plat=2*N_var[cutoff_index],
            # var_losecorr=N_var_losecorr[box_size_index]
        )
        assert np.isfinite(plateau)
        assert np.isfinite(plateau_std)

        # print(plateau.shape, plateau)
        plateau *= plateau_adjust
        plateau_std *= plateau_adjust
        # print(plateau.shape, plateau)
        assert np.isscalar(plateau_offset)
        plateau += plateau_offset
        # print(plateau.shape, plateau)

        # T_integrand     = T_integrand_func(N2_mean[box_size_index, :], plateau)
        # T_integrand_min = T_integrand_func(N2_mean[box_size_index, :] + N2_std[box_size_index, :], plateau, plateau_std)
        # T_integrand_max = T_integrand_func(N2_mean[box_size_index, :] - N2_std[box_size_index, :], plateau, plateau_std)
        T_integrand, T_integrand_unc = T_integrand_func(N2_mean[box_size_index, :], N2_std[box_size_index, :], plateau, plateau_std)
        T_integrand_max = T_integrand + T_integrand_unc
        T_integrand_min = T_integrand - T_integrand_unc
        # it's possible some of T_integrand_min are negative now so we prevent that
        T_integrand_min[T_integrand_min < 0] = 0

        # these are for the simple fit to C_N
        C_N_over_plateau = 1 - N2_mean[box_size_index, :] / (plateau) # countoscope paper eq. 1
        C_N_over_plateau_unc_sq = (N2_std[box_size_index, :] / (plateau))**2 + (N2_mean[box_size_index, :] * plateau_std / (2*plateau))**2
        C_N_over_plateau_unc = np.sqrt(C_N_over_plateau_unc_sq)

        M1_index, M2_index, thresh_line = get_M1_M2(file, L, t, T_integrand, plateau_source)
        print(thresh_line, show_thresh_line, ax)
        if thresh_line is not None and show_thresh_line and ax:
            ax.plot(t, thresh_line)
        else:
            if thresh_line is None:
                pass
                # assert False, f'you probably need to enter M1 and M2 for {file}'

        # full method
        (early_plot_x, early_plot_y), plot_late_integral_fit, D_of_L, D_of_L_min, D_of_L_max = timescaleintegral_full(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index)
        
        print(f'  final uncs {D_of_L_min/D_of_L:.3f} {D_of_L_max/D_of_L:.3f}')
        D_of_Ls[box_size_index] = D_of_L
        D_of_L_uncs[:, box_size_index] = [D_of_L - D_of_L_min, D_of_L_max - D_of_L]

        # fixed exponent
        (late_plot_x, late_plot_y), plot_late_integral_fit, D_of_L_fixexponent, D_of_L_min_fixexponent, D_of_L_max_fixexponent = timescaleintegral_fixexponent(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index)
        D_of_Ls_fixexponent[box_size_index] = D_of_L_fixexponent
        D_of_L_uncs_fixexponent[:, box_size_index] = [D_of_L_fixexponent - D_of_L_min_fixexponent, D_of_L_max_fixexponent - D_of_L_fixexponent]
        
        # nofit method
        # data_integral_full, data_integral_full_min, data_integral_full_max = do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, -1)
        
        # D_of_L_nofit, D_of_L_min_nofit, D_of_L_max_nofit = get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral_full, data_integral_full_min, data_integral_full_max, 0, 0, 0)
        # D_of_Ls_nofit[box_size_index] = D_of_L_nofit
        # D_of_L_uncs_nofit[:, box_size_index] = [D_of_L_nofit - D_of_L_min_nofit, D_of_L_max_nofit - D_of_L_nofit]
        
        # D_of_L_nofit, D_of_L_min_nofit, D_of_L_max_nofit = get_D_from_integral_contributions(L, 0, 0, 0, data_integral_full, data_integral_full_min, data_integral_full_max, 0, 0, 0)
        # D_of_Ls_nofit[box_size_index] = D_of_L_nofit
        # D_of_L_uncs_nofit[:, box_size_index] = [D_of_L_nofit - D_of_L_min_nofit, D_of_L_max_nofit - D_of_L_nofit]
        
        # nofit2 method
        D_of_L_nofit2, D_of_L_min_nofit2, D_of_L_max_nofit2 = timescaleintegral_nofit(t, L, T_integrand, T_integrand_min, T_integrand_max)
        D_of_Ls_nofit2[box_size_index] = D_of_L_nofit2
        D_of_L_uncs_nofit2[:, box_size_index] = [D_of_L_nofit2 - D_of_L_min_nofit2, D_of_L_max_nofit2 - D_of_L_nofit2]
        
        # simple fit
        D_simplefit, D_simplefit_unc, simplefit_func = C_N_simplefit(t, C_N_over_plateau, C_N_over_plateau_unc, L)
        D_of_Ls_simplefit[box_size_index] = D_simplefit
        D_of_L_uncs_simplefit[:, box_size_index] = D_simplefit_unc

        if box_size_indices:
            # the user chooses which boxes to display
            display = box_size_index in box_size_indices
        else:
            # we choose
            if num_boxes <= max_num_boxes:
                display = True
            else:
                divisor = math.ceil(num_boxes / max_num_boxes)
                display = box_size_index % divisor == 0

        if not ax:
            display = False

        if display:
            # N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L**2)) # countoscope eq. 2, countoscope overleaf doc
            # t_C_N_theory = np.logspace(-2, 5, 200)

            if rescale_x == RESCALE_X_L2:
                rescale_x_value = L**2
            else:
                rescale_x_value = 1

            if plot_color:
                color = plot_color
            elif colormap:
                color = colormap(box_size_index/len(box_sizes))
            else:
                color = common.colormap(box_size_index, 0, len(box_sizes))
            
            sep = sep_sizes[box_size_index]
            if sigma := data.get('particle_diameter'):
                label = f'$L={L/sigma:.2g}\sigma$'
                if sep_in_label:
                    label += f', $s={sep/sigma:.2g}\sigma$'
            else:
                label = f'L={L:.2f}'
                if sep_in_label:
                    label += f', s={sep:.2f}'
            if labels_on_plot: label = None

            # plot actual data
            if DONT_PLOT_ALL_POINTS_TO_REDUCE_FILESIZE and T_integrand.size > 100:
                points_to_plot = common.exponential_indices(t, num=1000)
            else:
                points_to_plot = np.index_exp[1:] # this is basically a functional way of writing points_to_plot = [1:]
                hide_data_after_nofit_cutoff = False
                print('i secretly disabled hide_data_after_nofit_cutoff') # without this we get an error because `points_to_plot[points_to_plot <= points_to_plot[above_nofit_cutoff]]` doesn't really make sense if points_to_plot is just an index_exp I think
                
            if hide_data_after_nofit_cutoff:
                above_nofit_cutoff = np.argmax(T_integrand[points_to_plot] < NOFIT_CROP_THRESHOLD/200)
                print('above nofit cutoff', above_nofit_cutoff)
                if above_nofit_cutoff: # it will return zero if none are above
                    print(points_to_plot)
                    print(above_nofit_cutoff)
                    # print(points_to_plot.shape, above_nofit_cutoff.shape)
                    points_to_plot = points_to_plot[points_to_plot <= points_to_plot[above_nofit_cutoff]]
            
            print(points_to_plot)
            plot_solid = points_to_plot[points_to_plot <= M2_index]
            plot_alpha = points_to_plot[points_to_plot >  M2_index]

            if plot_C_N_squared:
                ax.plot(t[plot_solid]/rescale_x_value, T_integrand[plot_solid], color=color, zorder=5, linestyle='none', marker='o', markersize=2, label=label)
                ax.plot(t[plot_alpha]/rescale_x_value, T_integrand[plot_alpha],  color=color, zorder=4, linestyle='none', marker='o', markersize=2, alpha=late_C_N_alpha)
            else:
                # assert False
                ax.plot(t[points_to_plot]/rescale_x_value, C_N_over_plateau[points_to_plot], color=color, zorder=5, linestyle='none', marker='o', markersize=2, label=label)
                
            # integrand_rescaled_ax.plot(t[1:]/L**2, T_integrand[1:])
            
            if show_short_fits:
                ax.plot(early_plot_x/rescale_x_value, early_plot_y, color=color, linestyle=':')
            if show_long_fits and plot_late_integral_fit:
                # ax.plot(t_fit_late/rescale_x_value, late_fit_func_real(t_fit_late, *late_popt), color=color, linestyle=':', linewidth=4)
                ax.plot(late_plot_x/rescale_x_value, late_plot_y, color=color, linestyle=':', linewidth=4)
            # else:
            #     assert False
            
            # full sDFT (provided D0)
            t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]))
            full_theory_N2 = lambda t, D0 : countoscope_theory.nmsd.inter_2d(t, D0, N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
            #                                                                some of these bits actually do nothing cause we take the ratio later
            full_theory_integrand = lambda t, D0 : (1 - 0.5 * full_theory_N2(t, D0) / (full_theory_N2(np.inf, D0)/2) )**2
            #                                                         ^^^^^^^^^^^^^^^^^^^^^^^^ should be N_var
            integrand_fit = full_theory_integrand(t_theory, D0)

            if show_theory:
                ax.plot(t_theory/rescale_x_value, integrand_fit, color='black', linewidth=1, label='sDFT inter. (no fit)' if box_size_index==0 else None)

            # this plots a straight line where we'll put the labels
            tt = t/rescale_x_value
            c = 3e-3
            m = 0.5
            label_line = c*tt**m
            # ax.plot(tt, label_line)

            # labels on plot
            if labels_on_plot:
                # t_index_for_text = int(t.size / 1000 * box_size_index**2)
                # t_index_for_text = int(M1_index/5)
                # t_index_for_text = int(10*L**1.5)
                # t_index_for_text = np.argmax(inte)
                # t_index_for_text = int(box_size_index**2.8 + 10)

                datavalues = T_integrand if plot_C_N_squared else C_N_over_plateau

                t_index_for_text = np.argmax(datavalues < label_line)

                angle = np.arctan(np.gradient(datavalues, t))[t_index_for_text] * 180/np.pi
                assert not np.isnan(angle)
                # plt.scatter(t_theory[t_index_for_text], N2_theory_points[t_index_for_text])
                
                if sigma and not np.isnan(sigma):
                    L_label = rf'$L={L/sigma:.2g}\sigma$'
                else:
                    L_label = rf'$L={L:.2g}$'

                used_color = labels_on_plot_font_color if labels_on_plot_font_color else color
                
                ax.text(
                    t[t_index_for_text]*LABELS_ON_PLOT_X_SHIFT, datavalues[t_index_for_text]*LABELS_ON_PLOT_Y_SHIFT, L_label,
                    horizontalalignment='center', color=used_color, fontsize=9,
                    transform_rotates_text=True, rotation=angle, rotation_mode='anchor',
                    zorder=10,
                )

            if show_T_of_L:
                T_of_L = countoscope_theory.timescaleint.T_of_L([L], D0, phi, sigma)
                ax.vlines(T_of_L, *ax.get_ylim(), color=color)

    if ax:
        if show_slope:
            slope = -1
            x1, y1 = 1.5e0, 3e-2
            slope_line_size = 12
            x2, y2 = x1*slope_line_size, y1*slope_line_size**slope
            ax.plot([x1, x2], [y1, y2], color=common.FIT_COLOR)
            ax.text(x1, y2, f'$t^{{{slope}}}$', ha='left', va='center', color=common.FIT_COLOR)
                    
        if show_legend:
            ax.legend(fontsize=legend_fontsize, loc='lower left')
        ax.loglog()
        if rescale_x == RESCALE_X_L2:
            ax.set_xlabel('$t/L^2$')
        else:
            ax.set_xlabel('$t$ (s)')

        if plot_C_N_squared:
            if not disable_ylabel: ax.set_ylabel(r'$C_N(t)^2$')
            ax.set_ylim(1e-5, 1.2e0)
        else:
            if not disable_ylabel: ax.set_ylabel(r'$C_N(t)$')
            ax.set_ylim(5e-4, 1.2e0)
        
        if show_nofit_cutoff:
            ax.hlines(NOFIT_CROP_THRESHOLD, *ax.get_xlim(), label='nofit crop threshold', color='gray')

    if save_data:
        common.save_data(f'visualisation/data/Ds_from_timescaleint_{plateau_source}{plateau_source_suffix}_{file}',
                Ds=D_of_Ls, D_uncs=D_of_L_uncs, Ls=box_sizes,
                particle_diameter=sigma, max_time_hours=data.get('max_time_hours'),)
        common.save_data(f'visualisation/data/Ds_from_timescaleint_fixexponent_{plateau_source}{plateau_source_suffix}_{file}',
                Ds=D_of_Ls_fixexponent, D_uncs=D_of_L_uncs_fixexponent, Ls=box_sizes,
                particle_diameter=sigma, max_time_hours=data.get('max_time_hours'),)
        common.save_data(f'visualisation/data/Ds_from_timescaleint_nofit_{plateau_source}{plateau_source_suffix}_{file}',
                Ds=D_of_Ls_nofit, D_uncs=D_of_L_uncs_nofit, Ls=box_sizes,
                particle_diameter=sigma, max_time_hours=data.get('max_time_hours'),)
        common.save_data(f'visualisation/data/Ds_from_timescaleint_nofit_cropped_{plateau_source}{plateau_source_suffix}_{file}',
                Ds=D_of_Ls_nofit2, D_uncs=D_of_L_uncs_nofit2, Ls=box_sizes,
                particle_diameter=sigma, max_time_hours=data.get('max_time_hours'),)
        # common.save_data(f'visualisation/data/Ds_from_timescaleint_nofit_cropped_noshort_{plateau_source}{plateau_source_suffix}_{file}',
        #         Ds=D_of_Ls_nofit_noshortfit, D_uncs=D_of_L_uncs_nofit_noshortfit, Ls=box_sizes,
        #         particle_diameter=sigma, max_time_hours=data.get('max_time_hours'),)
    common.save_data(f'visualisation/data/Ds_from_C_N_simplefit_{plateau_source}{plateau_source_suffix}_{file}',
                Ds=D_of_Ls_simplefit, D_uncs=D_of_L_uncs_simplefit, Ls=box_sizes,
                particle_diameter=sigma, max_time_hours=data.get('max_time_hours'),)

def timescaleintegral_nofit(t, L, T_integrand, T_integrand_min, T_integrand_max):
    early_fit_func, early_popt, early_integral, early_integral_min, early_integral_max = do_early_integral(f, t, T_integrand, T_integrand_min, T_integrand_max)

    stop_point = np.argmax(T_integrand < NOFIT_CROP_THRESHOLD)
    if stop_point == 0: # argmax returns zero if condition not met
        stop_point = -1
    data_integral_full2, data_integral_full2_min, data_integral_full2_max = do_data_integral(t[:stop_point], T_integrand[:stop_point], T_integrand_min[:stop_point], T_integrand_max[:stop_point], -1)

    D_of_L_nofit2, D_of_L_min_nofit2, D_of_L_max_nofit2 = get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral_full2, data_integral_full2_min, data_integral_full2_max, 0, 0, 0)
    return D_of_L_nofit2, D_of_L_min_nofit2, D_of_L_max_nofit2

def timescaleintegral_full(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index):
    early_fit_func, early_popt, early_integral, early_integral_min, early_integral_max = do_early_integral(f, t, T_integrand, T_integrand_min, T_integrand_max)
    
    early_for_plotting_t = np.logspace(-2, np.log10(t[1]))
    early_for_plotting_y = early_fit_func(early_for_plotting_t, *early_popt)
    
    data_integral, data_integral_min, data_integral_max = do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, M2_index)
    if L < 3:
            # the late integral fails here cause it's fully in the noise
            # but it doesn't matter for small boxes cause the noise is so low
        late_integral, late_integral_min, late_integral_max = 0, 0, 0
        plot_late_integral_fit = False
        
    else:
        try:
            late_integral, late_integral_min, late_integral_max, late_popt, late_fit_func_real, t_fit_late = do_late_integral(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index)
            
            print(f'  late integral contribution {late_integral/(early_integral+data_integral+late_integral):.3f}')
            plot_late_integral_fit = True

        except PositiveSlopeInTimescaleIntegralFitException:
            print('  full timescale integral failed')
            late_integral, late_integral_min, late_integral_max = 0, 0, 0
            plot_late_integral_fit = False

    D_of_L, D_of_L_min, D_of_L_max = get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral, data_integral_min, data_integral_max, late_integral, late_integral_min, late_integral_max)
    return (early_for_plotting_t, early_for_plotting_y), plot_late_integral_fit,D_of_L,D_of_L_min,D_of_L_max

def timescaleintegral_fixexponent(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index):
    early_fit_func, early_popt, early_integral, early_integral_min, early_integral_max = do_early_integral(f, t, T_integrand, T_integrand_min, T_integrand_max)
    
    data_integral, data_integral_min, data_integral_max = do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, M2_index)
    
    if L < 3:
            # the late integral fails here cause it's fully in the noise
            # but it doesn't matter for small boxes cause the noise is so low
        late_integral_fixexponent, late_integral_min_fixexponent, late_integral_max_fixexponent = 0, 0, 0
        plot_late_integral_fit = False
        late_plot_x, late_plot_y = None, None
        
    else:
        try:
            late_integral_fixexponent, late_integral_min_fixexponent, late_integral_max_fixexponent, late_popt_fixexponent, late_fit_func_real_fixexponent, late_plot_x = do_late_integral(t, L, T_integrand, T_integrand_min, T_integrand_max, M2_index, M1_index, fix_slope=True)
            late_plot_y = late_fit_func_real_fixexponent(late_plot_x, *late_popt_fixexponent)
            plot_late_integral_fit = True

        except PositiveSlopeInTimescaleIntegralFitException:
            print('  full timescale integral failed')
            plot_late_integral_fit = False
            late_plot_x, late_plot_y = None, None

    D_of_L, D_of_L_min, D_of_L_max = get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral, data_integral_min, data_integral_max, late_integral_fixexponent, late_integral_min_fixexponent, late_integral_max_fixexponent)
    return (late_plot_x, late_plot_y), plot_late_integral_fit, D_of_L, D_of_L_min, D_of_L_max
    

    # if title != None: # not if title b/c we want title='' to go down this path
    #     ax.set_title(title)
    # else:
    #     ax.set_title(fr'{file} timescale integrand, $\phi={phi:.3f}$, $\sigma={sigma}$, plateau:{var_label}', fontsize=10)


if __name__ == '__main__':  
    for file in common.files_from_argv('box_counting/data', 'counted_'):
        integ_fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        plateau_source = PLATEAU_SOURCE
        # plateau_source = 'sDFT'
        go(file, ax=ax, plateau_source=plateau_source, show_theory=SHOW_THEORY,
           labels_on_plot=LABELS_ON_PLOT, max_num_boxes=MAX_NUM_BOXES, show_short_fits=SHOW_TIMESCALEINTEGRAL_FIT_SHORT, show_long_fits=SHOW_TIMESCALEINTEGRAL_FIT_LONG,
           save_data=True, show_legend=SHOW_LEGEND, rescale_x=RESCALE_X, show_thresh_line=SHOW_THRESH_LINE,
           show_nofit_cutoff=SHOW_NOFIT_CUTOFF)#, plateau_adjust=1.05)
        
        # ax.set_xlim(5e-1, 1e5)

        if RESCALE_X == RESCALE_X_L2:
            ax.set_xlim(0.5e-1, 1e3)
            ax.set_ylim(1e-3, 1)

        ax.set_xlim(0.5, 1e6)
        ax.set_title(f'{file}, plateau: {plateau_source}')
        
        common.save_fig(integ_fig, f'box_counting/figures_png/integrand_{file}.png', dpi=300)