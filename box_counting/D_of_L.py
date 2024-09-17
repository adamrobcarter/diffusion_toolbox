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

SHOW_THEORY = False

METHOD_FULL = 0
METHOD_NOFIT = 1
# METHOD = METHOD_NOFIT
METHOD = METHOD_FULL

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
warnings.filterwarnings('ignore')

def do_early_integral(f, t, T_integrand, T_integrand_min, T_integrand_max):
    interval_end_index = 4
    early_fit_func = lambda t, a: f(t/a)**4
    early_fit_xdata = [0, t[1]]
    early_fit_ydata = [1, T_integrand[1]]

    early_popt, early_pcov = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata)
    early_integral = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt[0])[0]
        
    early_fit_ydata_min = [1, T_integrand_min[1]]
    early_popt_min, early_pcov_min = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata_min)
    early_integral_oneway = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt_min[0])[0]
        
    early_fit_ydata_max = [1, T_integrand_max[1]]
    early_popt_max, early_pcov_max = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata_max)
    early_integral_otherway = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt_max[0])[0]
        
    early_integral_min = min(early_integral, early_integral_oneway, early_integral_otherway)
    early_integral_max = max(early_integral, early_integral_oneway, early_integral_otherway)

    assert early_integral_min <= early_integral
    assert early_integral <= early_integral_max
    
    return early_fit_func,early_popt,early_integral,early_integral_min,early_integral_max

def do_late_integral(t, L, T_integrand, M2_index, M1_index):
    late_fit_func = lambda t, a, b: (a - b*t) - np.log(1 + np.exp(a - b*t))
    late_fit_xdata =     np.log(t          [M1_index:M2_index])
    late_fit_ydata = 2 * np.log(T_integrand[M1_index:M2_index])
    assert late_fit_xdata.size > 3
    late_popt, late_pcov = scipy.optimize.curve_fit(late_fit_func, late_fit_xdata, late_fit_ydata)
    a, b = late_popt
    a_unc = np.sqrt(late_pcov)[0, 0]
    b_unc = np.sqrt(late_pcov)[1, 1]
    late_fit_func_real = lambda t, a, b: np.exp(0.5 * late_fit_func(np.log(t), a, b))
    t_fit_late = np.logspace(np.log10(t[M1_index]), np.log10(t.max()))
    assert np.all(late_fit_func_real(t_fit_late, *late_popt) > 0), f'some of late_fit_func_real were negative for L={L}'
    late_integral = scipy.integrate.quad(late_fit_func_real, t[M2_index], t.max(), args=(a, b))[0]
            # late_integral_00 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t.max(), args=(a+a_unc))[0]
            # late_integral_01 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t.max(), args=tuple(late_popt))[0]
            # late_integral_10 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t.max(), args=tuple(late_popt))[0]
            # late_integral_11 = scipy.integrate.quad(late_fit_func_real, t[M2_index], t.max(), args=tuple(late_popt))[0]
            # print('isneg', late_fit_func_real(np.inf, *late_popt))
            # print('isneg', late_fit_func_real(1e10, *late_popt))
    t_to_inf = np.logspace(np.log10(t[M1_index]), 15, 1000)
    lf = late_fit_func_real(t_to_inf, *late_popt)
    assert np.all(lf >= 0)
    return late_integral, late_popt, late_fit_func_real, t_fit_late

def do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, M2_index):
    data_integral = scipy.integrate.trapezoid(T_integrand    [1:M2_index], t[1:M2_index])
    min_and_max = [
            scipy.integrate.trapezoid(T_integrand    [1:M2_index], t[1:M2_index]),
            scipy.integrate.trapezoid(T_integrand_max[1:M2_index], t[1:M2_index]),
            scipy.integrate.trapezoid(T_integrand_min[1:M2_index], t[1:M2_index])
        ]
    data_integral_min = min(min_and_max)
    data_integral_max = max(min_and_max)
    assert data_integral_min <= data_integral
    assert data_integral_max >= data_integral
    if data_integral_min == data_integral:
        warnings.warn('slightly strange uncertainty behavoir')
    return data_integral,data_integral_min,data_integral_max

def get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral, data_integral_min, data_integral_max, late_integral):
    total_integral     = early_integral     + data_integral     + late_integral
    total_integral_min = early_integral_min + data_integral_min + late_integral
    total_integral_max = early_integral_max + data_integral_max + late_integral
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

    assert D_of_L_min <= D_of_L
    assert D_of_L_max >= D_of_L
    return D_of_L, D_of_L_min, D_of_L_max


def C_N_simplefit(t, C_N_over_VarN, C_N_over_VarN_unc, L):
    
    func = lambda t, D : countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2

    fitting_points = common.exponential_integers(1, t.size//2)
    p0 = [0.05]
    popt, pcov = scipy.optimize.curve_fit(func, t[fitting_points], C_N_over_VarN[fitting_points], sigma=C_N_over_VarN_unc[fitting_points], p0=p0, maxfev=2000)
    
    D_from_fit = popt[0]
    D_from_fit_unc = np.sqrt(pcov[0][0])
    return D_from_fit, D_from_fit_unc, func

for file in common.files_from_argv('box_counting/data', 'counted_'):

    D_fig = plt.figure()
    T_fig = plt.figure()
    D_ax = D_fig.gca()
    T_ax = T_fig.gca()

    # integrand_fig = plt.figure()
    # integrand_ax = integrand_fig.gca()
    integrand_fig, (integrand_ax_old, mid_ax, extra_ax) = plt.subplots(3, 1, figsize=(6, 15))
    integ_fig, integ_axs = plt.subplots(1, 1, figsize=(6, 5))

    integrand_rescaled_fig = plt.figure()
    integrand_rescaled_ax = integrand_rescaled_fig.gca()

    C_N_fig, C_N_ax = plt.subplots(1, 1)

    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean   = data['N2_mean']
    N2_std    = data['N2_std']
    sigma     = data['particle_diameter']
    phi       = data['pack_frac']
    time_step = data['time_step']

    box_sizes = data['box_sizes']
    N_mean    = data['N_mean']
    N_var     = data['N_var']
    N_var_mod = data['N_var_mod']
    sep_sizes = data['sep_sizes']

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t = np.arange(0, num_timesteps) * time_step

    D_of_Ls     = np.full((num_boxes), np.nan)
    D_of_L_uncs = np.full((2, num_boxes), np.nan)
    T_of_Ls     = np.full((num_boxes), np.nan)
    
    D_of_Ls_nofit     = np.full((num_boxes), np.nan)
    D_of_L_uncs_nofit = np.full((2, num_boxes), np.nan)
    
    D_of_Ls_simplefit     = np.full((num_boxes), np.nan)
    D_of_L_uncs_simplefit = np.full((2, num_boxes), np.nan)


    for box_size_index in tqdm.trange(num_boxes, desc='timescale integral'):
        # T_of_Ls[box_size_index] = common.T_of_L(N2_mean[box_size_index, :], N_var[box_size_index], t)
        L = box_sizes[box_size_index]
        # D_of_Ls[box_size_index] = common.D_of_L(N2_mean[box_size_index, :], N_var[box_size_index], t, L) / D0


        print(f'L={L:.1f} N_var={N_var[box_size_index]:.4f}, N_var/L^2={N_var[box_size_index]/L**2:.4f}')


        plateau, plateau_std = box_counting.msd_single.get_plateau(N2_mean[box_size_index, :], file)

        # var_label = 'original'
        # var = N_var[box_size_index]
        var_label = 'plateau'
        var = plateau / 2

        T_integrand_func = lambda nmsd: (1 - 0.5 * nmsd / var )**2 # countoscope paper eq. 8
        T_integrand     = T_integrand_func(N2_mean[box_size_index, :])
        T_integrand_min = T_integrand_func(N2_mean[box_size_index, :] + N2_std[box_size_index, :])
        T_integrand_max = T_integrand_func(N2_mean[box_size_index, :] - N2_std[box_size_index, :])

        C_N_over_VarN = 1 - N2_mean[box_size_index, :] / (plateau) # countoscope paper eq. 1
        C_N_over_VarN_unc_sq = (N2_std[box_size_index, :] / (plateau))**2 + (N2_mean[box_size_index, :] * plateau_std / (2*plateau))**2
        C_N_over_VarN_unc = np.sqrt(C_N_over_VarN_unc_sq)

        # data integral setup
        if METHOD == METHOD_FULL:
            MIN_M2_INDEX = 8
        
            M2_index = None
            M2_index = None
            M1_index = None
            D0 = None

            if file.startswith('eleanorlong') or file.startswith('eleanor0.01') or file.startswith('brennan'):
                M2_index = min(int(300 * L), t.shape[0]-1)
                M2_index = max(M2_index, MIN_M2_INDEX)
                M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
                if file.startswith('eleanorlong'):
                    D0 = 0.028
                elif file.startswith('eleanor0.01') or file.startswith('brennan'):
                    D0 = 0.0416
            elif file == 'alice0.02':
                M2_index = min(int(60 * L), t.shape[0]-1)
                M2_index = max(M2_index, MIN_M2_INDEX)
                M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
                D0 = 0.0416
            elif file == 'alice0.02_overlapped3':
                M2_index = min(int(100 * L), t.shape[0]-1)
                M2_index = max(M2_index, MIN_M2_INDEX)
                M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
                D0 = 0.0416
            elif file == 'alice0.34':
                M2_index = min(int(40 * L), t.shape[0]-1)
                M2_index = max(M2_index, MIN_M2_INDEX)
                M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
                D0 = 0.0310
            elif file == 'alice0.66':
                M2_index = min(int(200 * L), t.shape[0]-1)
                M2_index = max(M2_index, MIN_M2_INDEX)
                M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
                D0 = 0.0175
            else:
                raise Exception('you need to define the parameters for this dataset')

            assert M2_index - M1_index > 3, f'M2_index = {M2_index}, M1_index = {M1_index}, t.shape[0]-1 = {t.shape[0]-1}'

            # ^^^^ this is the fundamental problem
        elif METHOD == METHOD_NOFIT:
            M2_index = -1 # bit hacky


        # full method
        early_fit_func, early_popt, early_integral, early_integral_min, early_integral_max = do_early_integral(f, t, T_integrand, T_integrand_min, T_integrand_max)
        data_integral, data_integral_min, data_integral_max = do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, M2_index)
        late_integral, late_popt, late_fit_func_real, t_fit_late = do_late_integral(t, L, T_integrand, M2_index, M1_index)
        
        D_of_L, D_of_L_min, D_of_L_max = get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral, data_integral_min, data_integral_max, late_integral)
        D_of_Ls[box_size_index] = D_of_L
        D_of_L_uncs[:, box_size_index] = [D_of_L - D_of_L_min, D_of_L_max - D_of_L]
        
        # nofit method
        data_integral_full, data_integral_full_min, data_integral_full_max = do_data_integral(t, T_integrand, T_integrand_min, T_integrand_max, -1)

        D_of_L_nofit, D_of_L_min_nofit, D_of_L_max_nofit = get_D_from_integral_contributions(L, early_integral, early_integral_min, early_integral_max, data_integral_full, data_integral_full_min, data_integral_full_max, 0)
        D_of_Ls_nofit[box_size_index] = D_of_L_nofit
        D_of_L_uncs_nofit[:, box_size_index] = [D_of_L_nofit - D_of_L_min_nofit, D_of_L_max_nofit - D_of_L_nofit]

        # simple fit
        D_simplefit, D_simplefit_unc, simplefit_func = C_N_simplefit(t, C_N_over_VarN, C_N_over_VarN_unc, L)
        D_of_Ls_simplefit[box_size_index] = D_simplefit
        D_of_L_uncs_simplefit[:, box_size_index] = D_simplefit_unc

        # if box_size_index%1 == 0:
        max_num_boxes = 12

        if num_boxes <= max_num_boxes:
            display = True
        else:
            divisor = math.ceil(num_boxes / max_num_boxes)
            display = box_size_index % divisor == 0


        if display:

            integrand_ax = integ_axs

            # print(f'early data late {early_integral:.1f}, {data_integral:.1f}, {late_integral:.1f}')
            # print('ratio', late_integral/late_integral1)
            N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L**2)) # countoscope eq. 2, countoscope overleaf doc
            t_C_N_theory = np.logspace(-2, 5, 200)
            
            integrand_ax.set_xlim(0.1, 1e5)


            color = common.colormap(box_size_index, 0, len(box_sizes))
            
            sep = sep_sizes[box_size_index]
            label = f'L={L:.2f}, s={sep:.2f}'
            if sigma := data.get('particle_diameter'):
                label = f'$L/\sigma={L/sigma:.2g}$, $s/\sigma={sep/sigma:.2g}$'
            line = integrand_ax.plot(t[1:M2_index], T_integrand[1:M2_index], color=color, label=label, zorder=5)
            integrand_ax.plot(t[M2_index:], T_integrand[M2_index:], alpha=0.5, color=color, zorder=4)
            integrand_rescaled_ax.plot(t[1:]/L**2, T_integrand[1:])
            
            t_fit_early = np.logspace(-2, np.log10(t[1]))
            integrand_ax.plot(t_fit_early, early_fit_func(t_fit_early, *early_popt), color=line[0].get_color(), linestyle=':')
            if METHOD == METHOD_FULL:
                integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, *late_popt), color=line[0].get_color(), linestyle=':', linewidth=4)
            
            
            integrand_ax.legend(fontsize=8)
            integrand_ax.loglog()
            integrand_ax.set_xlabel('$t$')
            

            # full sDFT (provided D0)
            if SHOW_THEORY:
                t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]))
                full_theory_N2 = lambda t, D0 : countoscope_theory.nmsd.inter_2d(t, D0, N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
                #                                                                some of these bits actually do nothing cause we take the ratio later
                full_theory_integrand = lambda t, D0 : (1 - 0.5 * full_theory_N2(t, D0) / (full_theory_N2(np.inf, D0)/2) )**2
                #                                                         ^^^^^^^^^^^^^^^^^^^^^^^^ should be N_var
                integrand_ax.plot(t_theory, full_theory_integrand(t_theory, D0), color='black', linewidth=1, label='sDFT inter. (no fit)' if box_size_index==0 else None)


            # full sDFT (fit to D0)
            # full_theory_N2 = lambda t, D0 : sDFT_interactions.sDFT_interactions(L=L, t=t, phi=phi, D0=D0, sigma=sigma)
            # full_theory_integrand = lambda t, D0 : (1 - 0.5 * full_theory_N2(t, D0) / (full_theory_N2(np.inf, D0)/2) )**2
            # #                                                         ^^^^^^^^^^^^^^^^^^^^^^^^ should be N_var
            # fitting_points = common.exponential_integers(1, M2_index)
            # popt_full, pcov_full = scipy.optimize.curve_fit(full_theory_integrand, t[fitting_points], T_integrand[fitting_points])
            # integrand_ax.plot(t_theory, full_theory_integrand(t_theory, *popt_full), color='black', linestyle='dotted', linewidth=1, label='sDFT inter. fit D' if box_size_index==0 else None)



            integrand_ax.set_ylim(1e-6, 1.2e0)
            # integrand_ax.set_ylim(1e-7, 1.2e0)

            C_N_ax.scatter(t, C_N_over_VarN, s=2)
            C_N_ax.plot(t, simplefit_func(t, D_simplefit), color='grey')


    D_ax.plot(box_sizes/sigma, D_of_Ls,  label=f'')
    T_ax.plot(box_sizes/sigma, T_of_Ls, label=f'')
    
    common.save_data(f'visualisation/data/Ds_from_timescaleint_{file}',
             Ds=D_of_Ls, D_uncs=np.zeros_like(D_of_L_uncs), Ls=box_sizes,
             particle_diameter=sigma)
    common.save_data(f'visualisation/data/Ds_from_timescaleint_nofit_{file}',
             Ds=D_of_Ls_nofit, D_uncs=np.zeros_like(D_of_L_uncs_nofit), Ls=box_sizes,
             particle_diameter=sigma)
    common.save_data(f'visualisation/data/Ds_from_C_N_simplefit_{file}',
             Ds=D_of_Ls_nofit, D_uncs=np.zeros_like(D_of_L_uncs_nofit), Ls=box_sizes,
             particle_diameter=sigma)
    common.save_data(f'visualisation/data/Ds_from_C_N_simplefit_{file}',
             Ds=D_of_Ls_simplefit, D_uncs=np.zeros_like(D_of_L_uncs_simplefit), Ls=box_sizes,
             particle_diameter=sigma)
    

    integrand_rescaled_ax.legend()
    integrand_rescaled_ax.loglog()
    integrand_rescaled_ax.set_xlabel('$t/L^2$')
    integrand_rescaled_ax.set_ylabel(r'integrand')
    integrand_rescaled_ax.set_ylim(1e-6, 1.1e0)

    # common.save_fig(integrand_rescaled_fig, f'box_counting/figures_png/integrand_rescaled_{file}.png', dpi=300)
    # integrand
    # integrand_ax.seintegrand_xlabel('$L/\sigma$')
    # integrand_ax.seintegrand_ylabel(r'$integrand(L)$')
    # integrand_ax.set_ylim(1e-4, 1.1e0)
    # integrand_fig.savefig('figures_png/integrand.png', dpi=300)

    integrand_ax.set_title(fr'{file} timescale integrand, $\phi={phi:.3f}$, $\sigma={sigma}$, var:{var_label}')
    common.save_fig(integ_fig, f'box_counting/figures_png/integrand_{file}.png', dpi=300)
    
    C_N_ax.semilogx()
    C_N_ax.semilogy()
    C_N_ax.set_ylim(1e-4, 1.1)
    common.save_fig(C_N_fig, f'box_counting/figures_png/C_N_{file}.png', dpi=300)