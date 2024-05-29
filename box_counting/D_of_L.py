import countoscope_theory.nmsd
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate, scipy.special, scipy.optimize, scipy.signal
import countoscope_theory.nmsd, countoscope_theory.structure_factor
import tqdm
import sys
import matplotlib.cm

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

for file in sys.argv[1:]:

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

    # D0 = { # countoscope paper, table 1
    #     0.02: 0.0416,
    #     0.34: 0.0310,
    #     0.66: 0.0175
    # }[phi]

    # data = common.load(f'box_counting/data/counted_dense_{file}.npz')
    data = common.load(f'box_counting/data/counted_{file}.npz')
    N2_mean = data['N2_mean']
    N2_std  = data['N2_std']
    # N_stats = data['N_stats']
    sigma = data['particle_diameter']
    phi = data['pack_frac']
    time_step    = data['time_step']

    # box_sizes = N_stats[:, 0]
    # N_mean    = N_stats[:, 1]
    # N_var     = N_stats[:, 2]
    box_sizes = data['box_sizes']
    N_mean    = data['N_mean']
    N_var     = data['N_var']
    sep_sizes = data['sep_sizes']

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t = np.arange(0, num_timesteps) * time_step

    D_of_Ls = np.full((num_boxes), np.nan)
    T_of_Ls = np.full((num_boxes), np.nan)

    # sDFT_data = common.load(f'data/sDFT_{phi}.npz')
    # late_integrals    = sDFT_data['late_integrals']
    # integrand_theorys = sDFT_data['integrand_theorys']
    # assert np.all(sDFT_data['box_sizes'] == box_sizes)

    for box_size_index in tqdm.trange(num_boxes, desc='timescale integral'):
        # T_of_Ls[box_size_index] = common.T_of_L(N2_mean[box_size_index, :], N_var[box_size_index], t)
        L = box_sizes[box_size_index]
        # D_of_Ls[box_size_index] = common.D_of_L(N2_mean[box_size_index, :], N_var[box_size_index], t, L) / D0

        print(f'L={L:.1f} N_var={N_var[box_size_index]:.4f}, N_var/L^2={N_var[box_size_index]/L**2:.4f}')

        T_integrand_func = lambda nmsd: (1 - 0.5 * nmsd / N_var[box_size_index] )**2 # countoscope paper eq. 8
        T_integrand     = T_integrand_func(N2_mean[box_size_index, :])
        T_integrand_min = T_integrand_func(N2_mean[box_size_index, :] + N2_std[box_size_index, :])
        T_integrand_max = T_integrand_func(N2_mean[box_size_index, :] - N2_std[box_size_index, :])
        # print(T_integrand_min.m

        C_N = N_var[box_size_index] - 0.5 * N2_mean[box_size_index, :] # countoscope paper eq. 1
        interval_end = (C_N[1] / N_var[box_size_index])**2
        interval_end_2 = T_integrand[1]
        # print(t[1], interval_end, interval_end_2)
        
        # early time integral
        interval_end_index = 4
        early_fit_func = lambda t, a: f(t/a)**4
        early_fit_xdata = [0, t[1]]
        early_fit_ydata = [1, T_integrand[1]]
        early_popt, early_pcov = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata)
        early_integral = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt[0])[0]
        
        early_fit_ydata = [1, T_integrand_min[1]]
        early_popt, early_pcov = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata)
        early_integral_min = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt[0])[0]
        
        early_fit_ydata = [1, T_integrand_max[1]]
        early_popt, early_pcov = scipy.optimize.curve_fit(early_fit_func, early_fit_xdata, early_fit_ydata)
        early_integral_max = scipy.integrate.quad(early_fit_func, 0, t[1], args=early_popt[0])[0]
        
        # late time integral
        
        # log_derivative = np.gradient(T_integrand) / T_integrand
        # log_derivative = scipy.signal.medfilt(log_derivative, 499)
        # minima_indices, props = scipy.signal.find_peaks(-log_derivative)
        
        # late_fit_func = lambda t, a, b: (a - b*t) - np.log(1 + np.exp(a - b*t))
        # late_fit_func = lambda t, a, b: (a - 3.35*t) - np.log(1 + np.exp(a - 3.35*t))
        # ^^^^ this is not okay!!!!!!!!

        MIN_M2_INDEX = 8
        
        M2_index = None
        M2_index = None
        M1_index = None
        D0 = None

        if file.startswith('eleanorlong'):
            M2_index = min(int(300 * L), t.shape[0]-1)
            M2_index = max(M2_index, MIN_M2_INDEX)
            M1_index = M2_index // 2 # t is linear so we can just halve the index to halve the time
            D0 = 0.028
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

        # late_fit_xdata =     np.log(t          [M1_index:M2_index])
        # late_fit_ydata = 2 * np.log(T_integrand[M1_index:M2_index])
        # late_popt, late_pcov = scipy.optimize.curve_fit(late_fit_func, late_fit_xdata, late_fit_ydata)
        # late_fit_func_real = lambda t, a, b: np.exp(0.5 * late_fit_func(np.log(t), a, b))
        # t_fit_late = np.logspace(np.log10(t[M1_index]), np.log10(t.max()))
        # # t_fit_late = np.logspace(-2, np.log10(t.max()))
        # assert np.all(late_fit_func_real(t_fit_late, *late_popt) > 0), f'some of late_fit_func_real were negative for L={L}'
        # # late_integral = scipy.integrate.quad(late_fit_func_real, t[M2_index], t.max(), args=tuple(late_popt))[0]
        # # print('isneg', late_fit_func_real(np.inf, *late_popt))
        # # print('isneg', late_fit_func_real(1e10, *late_popt))
        # t_to_inf = np.logspace(np.log10(t[M1_index]), 15, 1000)
        # lf = late_fit_func_real(t_to_inf, *late_popt)
        # assert np.all(lf >= 0)

        # def diy_quad(x, y):
        #     assert x.shape == y.shape
        #     I = 0
        #     for i in range(len(x)-1):
        #         I += (x[i+1] - x[i]) * (y[i+1] + y[i]) / 2
        #     return I

        # late_integral_res = scipy.integrate.quad(late_fit_func_real, t[M2_index], 1e7, full_output=True, args=tuple(late_popt))
        # late_integral1 = late_integral_res[0]

        # TODO: CAN WE NOT JUST TAKE THE INTEGRAL ANALYTICALLY!!!!


        # late_integral = diy_quad(t_to_inf, lf)
        # extra_decade_t = np.logspace(np.log10(t_to_inf.max()), np.log10(t_to_inf.max())*10)
        # extra_decade_lf = late_fit_func_real(extra_decade_t, *late_popt)
        # extra_decade_int = diy_quad(extra_decade_t, extra_decade_lf)
        # assert extra_decade_int / late_integral < 0.01, f'extra {extra_decade_int}, int {late_integral}'

        # print(late_integral_res[2])
        # print(late_integral_res[1])

        # print(diy_quad(t_to_inf, lf))
        # print()

        
        # late_integral_res = scipy.integrate.quad(late_fit_func, np.log(t[M2_index]), np.log(1e7), full_output=True, args=tuple(late_popt))
        # late_integral = late_integral_res[0]
        # late_integral = np.exp(0.5 * late_integral)

        # print('ratio', late_integral/late_integral1)

        # assert late_integral >= 0, f'late_integral = {late_integral} for L={L}'

        # data integral
        data_integral     = scipy.integrate.trapezoid(T_integrand    [1:M2_index], t[1:M2_index])
        data_integral_min = scipy.integrate.trapezoid(T_integrand_min[1:M2_index], t[1:M2_index])
        data_integral_max = scipy.integrate.trapezoid(T_integrand_max[1:M2_index], t[1:M2_index])

        # scipy.integrate.quad(early_func, 0, t[1])
        # T_of_Ls[box_size_index] = 2 * (scipy.integrate.quad(fit_func, 0, t[1])[0] + scipy.integrate.trapezoid(T_int[1:], t[1:])) # countoscope eq. 8
        
    
        L_2 = L
        # N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc
        # integrand_theory = lambda t: (1 - 0.5 * N2_func_full(t, D0) / N_mean[box_size_index])**2
        # #                                                             ^^^^^^ this should be N_var
        # integrand_offset = integrand_theory(t[M2_index]) - T_integrand[M2_index]
        # integrand_offset = 0


        # N2_theory_inter = lambda t: 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t, phi, D0, sigma)
        # integrand_theory_inter = lambda t: (1 - 0.5 * N2_theory_inter(t) / (N2_theory_inter([np.inf])/2))**2
        #                                                                   ^^^^^^ this should be N_var
        # late_integral = scipy.integrate.quad(integrand_theory_inter, t[M2_index], np.inf)[0]
        # late_integral = late_integrals[box_size_index]
        # assert late_integral >= 0, f'late_integral = {late_integral} for L={L}'

        # late integral: countoscope SI fit
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
        # t_fit_late = np.logspace(-2, np.log10(t.max()))
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

        print(early_integral_min, early_integral, early_integral_max)
        print(data_integral_min, data_integral, data_integral_max)

        
        T_of_Ls[box_size_index] = 2 * (early_integral + data_integral + late_integral) # countoscope eq. 8
        # T_of_Ls[box_size_index] = 2 * (early_integral + data_integral) # countoscope eq. 8
        alpha_T = 0.561244
        D_of_Ls[box_size_index] = alpha_T * L**2 / (4 * T_of_Ls[box_size_index]) # countoscope eq. 9


        # popt_full, pcov_full = scipy.optimize.curve_fit(integrand_full, t[:M2_index], T_integrand[:M2_index], maxfev=2000)

    
        # if box_size_index%1 == 0:
        if True:

            integrand_ax = integ_axs

            # integrand_ax.set_title(rf'$\phi={phi}$')

            # N_var_theory = N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L=L, t=np.array([np.inf]), phi=phi, D0=D0, sigma=sigma)[0]
            # print('var', N_var_theory, N_var[box_size_index], N_var_theory / N_var[box_size_index])
            #              ^^^^ should be density * box volume?

            # print(f'early data late {early_integral:.1f}, {data_integral:.1f}, {late_integral:.1f}')
            # print('ratio', late_integral/late_integral1)
            N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L**2)) # countoscope eq. 2, countoscope overleaf doc
            t_C_N_theory = np.logspace(-2, 5, 200)
            # N2_theory = N2_func_full(t_C_N_theory, D0)
            # integrand_theory = (1 - 0.5 * N2_theory / N_var_theory)**2
            # integrand_theory2 = (1 - 0.5 * N2_theory / N_mean[box_size_index])**2
            # integrand_ax.plot(t_C_N_theory, integrand_theory, color='grey', linewidth=0.8, zorder=5)
            # integrand_ax.plot(t_C_N_theory, integrand_theory2, color='black', linewidth=0.5, zorder=5)
            # integrand_rescaled_ax.plot(t_C_N_theory/L**2, integrand_theory, color='black', linewidth=0.5, zorder=5, alpha=1)
            integrand_ax.set_xlim(0.1, 1e5)

            # integrand_ax.plot(t_C_N_theory, integrand_theory_inter(t_C_N_theory), color='black', linewidth=0.5, zorder=5, alpha=1)
            # integrand_ax.plot(t_C_N_theory, integrand_full(t_C_N_theory, D0, N_mean), color='black', linewidth=0.5, zorder=5, alpha=1.0)
            # print('nvar', popt_full[0]/N_var[box_size_index])

            # integrand_ax.plot(t_C_N_theory, f(4*D0*t_C_N_theory/L**2)**4, color='black', linewidth=0.5, zorder=5, alpha=1.0, linestyle='dotted')
            
            

            #T_integrand = (1 - 0.5 * N2_mean[box_size_index, :] / N_var[box_size_index] )**2 # countoscope paper eq. 8
            # mid_ax.plot(t[1:], (N2_mean[box_size_index, :] / N_var[box_size_index])[1:], label=f'L={L:.1f}')
           
            # mid_ax.plot(t[1:],            1/2 * (N2_mean[box_size_index, 1:] / N_var[box_size_index]), label=f'L={L:.2f}')
            # mid_ax.plot(t_C_N_theory[1:], 1/2 * (N2_theory[1:]               / N_var[box_size_index]), color='black', linewidth=0.5)
            
            # mid_ax.plot(t[1:],            1 - 1/2 * (N2_mean[box_size_index, 1:] / N_mean[box_size_index]), label=f'L={L:.1f}')
            # mid_ax.plot(t_C_N_theory[1:], 1 - 1/2 * (N2_theory[1:]               / N_mean[box_size_index]), color='black', linewidth=0.5)
            
            # mid_ax.plot(t[1:],            ( 1 - 1/2 * (N2_mean[box_size_index, 1:] / N_mean[box_size_index]) )**2, label=f'L={L:.1f}')
            # mid_ax.plot(t_C_N_theory[1:], ( 1 - 1/2 * (N2_theory[1:]               / N_mean[box_size_index]) )**2, color='black', linewidth=0.5)
            
            # mid_ax.loglog()
            # print(f'L={L}, N_mean/N_var={N_mean[box_size_index]/N_var[box_size_index]}')


            color =  matplotlib.cm.afmhot(np.interp(box_size_index, (0, len(box_sizes)), (0.2, 0.75)))
            
            sep = sep_sizes[box_size_index]
            label = f'L={L:.2f}, s={sep:.2f}'
            line = integrand_ax.plot(t[1:M2_index], T_integrand[1:M2_index], color=color, label=label, zorder=5)
            integrand_ax.plot(t[M2_index:], T_integrand[M2_index:], alpha=0.5, color=color, zorder=4)
            integrand_rescaled_ax.plot(t[1:]/L**2, T_integrand[1:])
            
            t_fit_early = np.logspace(-2, np.log10(t[1]))
            integrand_ax.plot(t_fit_early, early_fit_func(t_fit_early, *early_popt), color=line[0].get_color(), linestyle=':')
            
            # integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, *late_popt), color=line[0].get_color(), linestyle=':', linewidth=3)
            
            t_to_inf = np.logspace(np.log10(t[M1_index]), 5)
            # integrand_ax.plot(t_to_inf, integrand_theorys[box_size_index, :], color=line[0].get_color(), linestyle=':', linewidth=3)

            # integrand_ax.plot(t_C_N_theory, integrand_theory(t_C_N_theory), color=line[0].get_color(), linewidth=0.8, zorder=5, alpha=1)
            
            
            integrand_ax.legend(fontsize=8)
            integrand_ax.loglog()
            integrand_ax.set_xlabel('$t$')
            # print(L, *late_popt)
            
            # log_derivative = np.gradient(T_integrand) / T_integrand
            # log_derivative = scipy.signal.medfilt(log_derivative, 499)
            # minima_indices, props = scipy.signal.find_peaks(-log_derivative)
            # # integrand_ax.vlines(t[minima_indices[0]], -0.02, 1,    color=line[0].get_color(), linestyle=':')
            # extra_ax    .vlines(t[minima_indices[0]], -0.02, 0.01, color=line[0].get_color(), linestyle=':')

            # extra_ax.plot(t[1:], log_derivative[1:])
            # extra_ax.semilogx()

            # full sDFT (provided D0)
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



            integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, a, b), color='grey', linewidth=1, label='counto. SI fit' if box_size_index==0 else None, zorder=6)
            integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, a+a_unc, b+b_unc), color='grey', linewidth=1.5, zorder=6)
            integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, a+a_unc, b-b_unc), color='grey', linewidth=2, zorder=6)
            integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, a-a_unc, b+b_unc), color='grey', linewidth=2.5, zorder=6)
            integrand_ax.plot(t_fit_late, late_fit_func_real(t_fit_late, a-a_unc, b-b_unc), color='grey', linewidth=3, zorder=6)
            


            integrand_ax.set_ylim(1e-5, 1.2e0)
            # integrand_ax.set_ylim(1e-7, 1.2e0)


    D_ax.plot(box_sizes/sigma, D_of_Ls/D0,  label=f'')
    T_ax.plot(box_sizes/sigma, T_of_Ls, label=f'')
    common.save_data(f'visualisation/data/Ds_from_timescaleint_{file}',
             Ds=D_of_Ls, D_uncs=np.zeros_like(D_of_Ls), Ls=box_sizes,
             particle_diameter=sigma)


    D_ax.legend()
    D_ax.loglog()
    D_ax.set_xlabel('$L/\sigma$')
    D_ax.set_ylabel(r'$D(L) / D_0$')
    common.save_fig(D_fig, f'box_counting/figures_png/D_of_L_{file}.png')

    T_ax.legend()
    T_ax.loglog()
    T_ax.set_xlabel('$L/\sigma$')
    T_ax.set_ylabel(r'$T(L)$')
    common.save_fig(T_fig, f'box_counting/figures_png/T_of_L_{file}.png')

    integrand_rescaled_ax.legend()
    integrand_rescaled_ax.loglog()
    integrand_rescaled_ax.set_xlabel('$t/L^2$')
    integrand_rescaled_ax.set_ylabel(r'integrand')
    integrand_rescaled_ax.set_ylim(1e-6, 1.1e0)
    common.save_fig(integrand_rescaled_fig, f'box_counting/figures_png/integrand_rescaled_{file}.png', dpi=300)
    # integrand
    # integrand_ax.seintegrand_xlabel('$L/\sigma$')
    # integrand_ax.seintegrand_ylabel(r'$integrand(L)$')
    # integrand_ax.set_ylim(1e-4, 1.1e0)
    # integrand_fig.savefig('figures_png/integrand.png', dpi=300)
    integrand_ax.set_title(fr'{file} timescale integrand, $\phi={phi:.3f}$, $\sigma={sigma}$')
    common.save_fig(integ_fig, f'box_counting/figures_png/integrand_{file}.png', dpi=300)

    # common.save_fig(integ_fig,  f'box_counting/figures_png/integrand_{file}.png', dpi=300)