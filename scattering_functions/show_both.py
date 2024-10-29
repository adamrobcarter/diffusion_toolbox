import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings, math
import scipy.optimize, scipy.signal
import matplotlib.cm
import visualisation.Ds_overlapped

subplot_i = 0

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)
# target_ks = list(np.logspace(np.log10(0.02), np.log10(8), 20))
# target_ks = list(np.logspace(np.log10(0.02), np.log10(0.1), 25))

SHOW_SEGRE_PUSEY_RESCALING_AXIS = False

D_ERROR_THRESH = 10 # point ignored if D_unc/D > this. was 0.1

EXP_FIT = 0
T_MINUS_T0_FIT = 1
EXP_TIMES_CONST_FIT = 2
DOMINIGUEZ_FIT = 3

FIT = EXP_TIMES_CONST_FIT
# FIT = DOMINIGUEZ_FIT

DO_SHORT_LINEAR_FIT = False
FIT_WITH_ZERO_POINT = True

CROP_AT_SK_MINIMUM = False

SHOW_TWOSTAGE_FIT = False

FILTER_RISING_MULTIPLE = 2.5

def show_single_F_type(
        file, type_index, Ftype, fig, axes,
        num_displayed_ks, mult=False, do_fits=True,
        markersize=5, errorbar_alpha=0.2, bad_alpha=0.2, color=None
    ):
    lin_short_axes, lin_axes, log_axes, D_axes = axes
    # lin_short_axes, lin_axes, log_axes, D_axes, extra_axes = axes

    Ds_for_saving = []
    D_uncs_for_saving = []
    ks_for_saving = []
        
    Ds_for_saving_short = []
    D_uncs_for_saving_short = []
    ks_for_saving_short = []

    Ds_for_saving_long = []
    D_uncs_for_saving_long = []
    ks_for_saving_long = []

    Ds_for_saving_first = []
    D_uncs_for_saving_first = []
    ks_for_saving_first = []

    Ds_for_saving_D1 = []
    D_uncs_for_saving_D1 = []
    ks_for_saving_D1 = []

    Ds_for_saving_D2 = []
    D_uncs_for_saving_D2 = []
    ks_for_saving_D2 = []

    if Ftype == 'f' and SHOW_SEGRE_PUSEY_RESCALING_AXIS:
        sp_fig, sp_axs = plt.subplots(2, 1, figsize=(4, 6))

    load = 'F' if Ftype=='f' else Ftype
    d = common.load(f"scattering_functions/data/{load}_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]
    particle_diameter = d['particle_diameter']

    ks_2D = len(k_all.shape) == 2

    num_ks = k_all.shape[1] if ks_2D else k_all.size

    every_nth_k = int(math.ceil(num_ks / num_displayed_ks))
    every_nth_k = max(every_nth_k, 1)

    graph_i = 0

    if CROP_AT_SK_MINIMUM:
        S_of_k = F_all[0, :]
        min_index = np.nanargmin(S_of_k)
        k_at_min = k_all[0, min_index]

    try:
        D0, _, _ = visualisation.Ds_overlapped.get_D0(file)
    except Exception as e:
        print(e)
        D0 = None

    def print_D(D):
        if D0 == None:
            return f'$D={D:.3g}\mathrm{{\mu m^2/s}}$'
        else:
            return f'$D={D/D0:.2g}D_0$'

    if color == None:
        color = colors[type_index]

    for k_index in range(num_ks):
        if Ftype == 'DDM' or (not ks_2D):
            ks = k_all
        else:
            ks = k_all[0, :]

        # k_index = np.argmax(ks >= target_k)

        k = ks[k_index]


        # print(f'k: target {target_k:.3f}, real {k:.3f}, index {k_index}, 2pi/k={2*np.pi/k:.1f}um')
        k_str = f'k {k:.3f}, index {k_index}, 2pi/k={2*np.pi/k:.1f}um'
        if particle_diameter:
            k_str += f'={2*np.pi/k/particle_diameter:.2g}σ'
        print(k_str)

        if CROP_AT_SK_MINIMUM and k < k_at_min:
            print('  skipping - below S(k) minimum')
            continue

        f     = F_all    [:, k_index]
        f_unc = F_unc_all[:, k_index]
        if Ftype == 'f':
            f /= F_all[0, k_index]
            f_unc_sq_all = (F_unc_all / F_all[0, :])**2 + (F_all * F_unc_all[0, :] / F_all[0, :]**2)**2
            # assert f_unc_sq.shape == f_unc.shape, f'{f_unc_sq.shape} {f_unc.shape}'
            f_unc = np.sqrt(f_unc_sq_all)[:, k_index]

            # assert k_index != 0, "k_index == 0, which probably means everything is 1 because of the normalisation"
            # if k_index == 0:
            #     warnings.warn("k_index == 0, which probably means everything is 1 because of the normalisation")
            #     continue

        if np.all(f == 0):
            print('  all f zero!')
        if np.all(f == 1):
            print('  all f == 1!')
        
        if np.isfinite(particle_diameter):
            k_label = fr'$k\sigma={k*particle_diameter:.2g}$'
            L_label = fr'$L\approx{2*np.pi/k:.3g}\mathrm{{\mu m}}={2*np.pi/k/particle_diameter:.2g}\sigma$'
        else:
            k_label = fr'$k={k:.2g}$'
            L_label = fr'$L\approx{2*np.pi/k:.3g}\mathrm{{\mu m}}$'

        label = fr"$k={k:.3f}\mathrm{{\mu m}}$ ({L_label})"

        if common.nanfrac(f) == 1:
            print(f'all nan at k={k:.2g} (i={graph_i})')
            continue

        
        display = False
        if k_index % every_nth_k == 0:
            display = True
            display = False
        if k_index < 10:
            display = True
        if display:

            ax = lin_axes[graph_i]
            ax_short = lin_short_axes[graph_i]
            D_ax = D_axes[graph_i]
            log_ax = log_axes[graph_i]
            log_ax.semilogy()
            # extra_ax = extra_axes[graph_i]

            graph_i += 1

        if do_fits:
            # first point D:
            if t.size == 1:
                first_index = 0
            else:
                assert t[0] == 0
                assert f[0] == 1
                first_index = 1
            if f[first_index] > 1e-2:
                D_first = -1 / (k**2 * t[first_index]) * np.log(f[first_index])
                D_unc_first = 1 / (k**2 * t[first_index] * f[first_index]) * f_unc[first_index]
                # print(f'  first D={D_first:.3g}')
                assert D_unc_first >= 0, f'D_unc_first={D_unc_first:.3f} = 1 / ({k:.3f}**2 * {t[1]:.3f} * {f[1]:.3f}) * {f_unc[1]:.3f}'
                Ds_for_saving_first    .append(D_first)
                D_uncs_for_saving_first.append(D_unc_first)
                ks_for_saving_first    .append(k)
            else:
                print('  skipped first')

        # if Ftype == 'f':
        #     noise_thresh = 1e-2 # for eleanorlong!!
        #     time_thresh = 200
        # elif Ftype == 'Fs':
        #     noise_thresh = 1.7e-2
        #     time_thresh = 400
        # elif Ftype == 'DDM':
        #     noise_thresh = 1e-3
        #     time_thresh = 400
        # f_noise   = f < noise_thresh
        # f_toolong = t > time_thresh
        #     # f_bad   = f_noise | f_toolong # f_toolong should depend on k!!!!
        # f_bad   = f_noise
        # f_bad[0] = True

        f_bad = np.full(f.shape, False)
        if t.size > 1:
            # new noise identification idea
            # first noise point is first point where the gradient is no longer getting more negative
            grad = np.gradient(np.log10(f), np.log10(t))
            width = None
            if file.startswith('marine'):
                prominence = 1
            elif file.startswith('eleanorlong'):
                prominence = 0.1
                if file == 'eleanorlong001':
                    prominence = 0.1
                if file == 'eleanorlong034':
                    prominence = 1
                if file == 'eleanorlong066':
                    prominence = 2
                # width = 1
            elif file.startswith('brennan'):
                prominence = 1 # 0.66
            else:
                prominence = 0.1 # was 0.01
            # prominance filters out peaks that are just noise. A smoothing filter would probably be better though
            # increase prominance to allow bigger noise
            peaks, props = scipy.signal.find_peaks(-grad, prominence=prominence, width=width)
            
            if len(peaks):
                f_bad[peaks[0]:] = True
            else:
                print('  no peaks in -grad found?!')
            

        # now we do another filter (should this be the only one?)
        # if the signal increases by more than `thresh` times, it's noise
        increase = f[1:] / f[:-1]
        above_thresh = increase > FILTER_RISING_MULTIPLE
        if above_thresh.sum():
            f_bad[np.argmax(above_thresh):] = True

        if np.sum(~f_bad) > 1:
            increase_good = f[~f_bad][1:] / f[~f_bad][:-1]
            print(f'  biggest good rise: {increase_good.max():.2f}')

        # negative points are definitely noise
        f_bad[f<=0] = True

        # don't use the zero point to fit
        if not FIT_WITH_ZERO_POINT:
            f_bad[0] = True
        else:
            assert f_bad[0] == False

        if display:
            
            # full time linear axis
            # good points, bad points
            ax.plot(t[~f_bad], f[~f_bad], color=color, linestyle='', label=f'{Ftype} {file}', marker='.', markersize=markersize)
            ax.plot(t[~f_bad], f[~f_bad], color=color, linestyle='', alpha=bad_alpha,               marker='.', markersize=markersize)
            # good errorbars, bad errorbars
            ax.errorbar(t[~f_bad], f[~f_bad], yerr=f_unc[~f_bad], color=color, linestyle='', alpha=errorbar_alpha,                 marker='none')
            ax.errorbar(t[ f_bad], f[ f_bad], yerr=f_unc[ f_bad], color=color, linestyle='', alpha=min(errorbar_alpha, bad_alpha), marker='none')

            # log log axis
            # good points, bad points
            log_ax.plot(t[~f_bad], f[~f_bad], color=color, linestyle='', label=f'{Ftype} {file}', marker='.', markersize=markersize)
            log_ax.plot(t[ f_bad], f[ f_bad], color=color, linestyle='', alpha=bad_alpha,               marker='.', markersize=markersize)
            # good errorbars, bad errobars
            log_ax.errorbar(t[~f_bad], f[~f_bad], yerr=f_unc[~f_bad], color=color, linestyle='', alpha=errorbar_alpha,                 marker='none')
            log_ax.errorbar(t[ f_bad], f[ f_bad], yerr=f_unc[ f_bad], color=color, linestyle='', alpha=min(errorbar_alpha, bad_alpha), marker='none')

            log_ax.set_ylim(f.min()/1.2, f.max()*1.2)

            # short time linear axis
            ax_short.plot(t, f, color=color, linestyle='', marker='.', markersize=markersize)
            ax_short.errorbar(t, f, yerr=f_unc, color=color, linestyle='', alpha=errorbar_alpha)
            
            SHORT_PLOT_MAX_X = 10 # for eleanor etc
            SHORT_PLOT_Y0 = f[0:int(SHORT_PLOT_MAX_X+1)].min() # for eleanor

            if file.startswith('marine'):
                SHORT_PLOT_MAX_X = 1 # for marine
                SHORT_PLOT_Y0 = 0 # for marine

            ax_short.set_xlim(0, SHORT_PLOT_MAX_X)
            ax_short.set_ylim(SHORT_PLOT_Y0, f[0:int(SHORT_PLOT_MAX_X+1)].max())


            
            # print(f'{Ftype} k={k:.2f} avg err = ')

                # ax.semilogx()
            ax.semilogy()
            newtitle = '\n' + Ftype + ':' + label
            if ax.get_title() != newtitle:
                ax.set_title(ax.get_title() + newtitle, fontsize=9)
                # print(~f_bad[1:5])
                # print(np.argmax(f_bad[1:5]))
                # print(np.argmax(f_bad[1:] ))
                # ax.set_xlim(0, max(t[np.argmax(f_bad[1:])], 100))
                # ax.set_xlim(0, 100)
                # ax.set_ylim(9.9e-1, 1.01)
                # ax.set_xlim(0, 1000)
                # ax.set_ylim(1e-3 * (1/k) , 1.1)
            end_plot_time = 1/(k+0.005) * 10 # WAS 200
            end_plot_y = f[np.argmax(t > end_plot_time)] * 0.99
            print('@@@', np.argmax(t > end_plot_time), f[np.argmax(t > end_plot_time)])
            if end_plot_y > 1: end_plot_y = 0.9
            if end_plot_y < 1e-4: end_plot_y = 1e-4
            ax.set_ylim(end_plot_y , 1.01)
                # ax.set_xlim(0, 1/k**2 * 100)
            ax.set_xlim(-end_plot_time/50, end_plot_time)

            # ax.set_ylim(min(max(func(t[np.argmax(f_bad[1:])], *f_popt), 1e-3), 9.9e-1), 1.01)
            # ax.set_ylim(1e-3, 1.1)



        negative = f <= 0
        if negative.sum() > 0:
            print(f'  negative: {negative.sum()/negative.size:.2f}')

        if do_fits:
            # fits
            # we fit to log(f(k, t)) because otherwise points near one make a much
            # larger impact to the fit than points near zero
            if FIT == T_MINUS_T0_FIT:
                func = lambda t, D, t0 : np.exp(-(t-t0) * k**2 * D)
                log_func = lambda t, D, t0: np.log10( func(t, D, t0) )
                # log_unc = lambda x, dx : 0.5 * np.log((x+dx)/(x-dx)) # what is this?
                # log_f_unc  = log_unc(f, f_unc )
            elif FIT == EXP_TIMES_CONST_FIT:
                func = lambda t, D, c : c * np.exp(-t * k**2 * D)
                log_func = lambda t, D, c: np.log10( func(t, D, c) )
                # log_unc = lambda x, dx : 0.5 * np.log((x+dx)/(x-dx)) # what is this?
                # log_f_unc  = log_unc(f, f_unc )
            elif FIT == EXP_FIT:
                func = lambda t, D : np.exp(-t * k**2 * D)
                log_func = lambda t, D: np.log10( func(t, D) )
                # log_unc = lambda x, dx : 0.5 * np.log((x+dx)/(x-dx)) # what is this?
                # log_f_unc  = log_unc(f, f_unc )
            elif FIT == DOMINIGUEZ_FIT:
                assert '034' in file
                Lh = 0.2*particle_diameter
                if 1/k < Lh:
                    func = lambda t, D : np.exp(-t * k**2 * D)
                    log_func = lambda t, D: np.log10( func(t, D) )
                else:
                    func = lambda t, D : np.exp(-t * k * D / Lh)
                    log_func = lambda t, D: np.log10( func(t, D) )
                
            p0 = [0.05]
            p0 = None
            

            if (~f_bad).sum() < 2:
                print('  no good data, skipping rest')
                continue
                
            # total fit
            f_popt, f_pcov = scipy.optimize.curve_fit(
                    log_func, t[~f_bad], np.log10(f[~f_bad]),
                    # p0=p0, sigma=log_f_unc[~f_bad], absolute_sigma=True
                )
            
            if FIT == EXP_TIMES_CONST_FIT:
                print(f'  total c = {f_popt[1]:.3g} pm {np.sqrt(f_pcov[1, 1])}')

            if FIT_WITH_ZERO_POINT:
                theory_start = 0.1
            else:
                theory_start = t[1]
            t_th = np.logspace(np.log10(theory_start), np.log10(t[-1]))

            if np.isfinite(np.sqrt(f_pcov)[0,0]):
                D = f_popt[0]
                D_unc = np.sqrt(f_pcov)[0][0]

                if display:
                    ax      .plot(t_th, func(t_th, *f_popt), color=color, linestyle='dashdot', label=f'total {print_D(D)}')
                    ax_short.plot(t_th, func(t_th, *f_popt), color=color, linestyle='dashdot', label=f'total {print_D(D)}')
                    log_ax  .plot(t_th, func(t_th, *f_popt), color=color, linestyle='dashdot', label=f'total {print_D(D)}')
                        # ax.plot(t_th, func(t_th, *Fs_popt), color='tab:orange',   linestyle='dotted')
                    D_ax.hlines(f_popt[0], t.min(), t.max(), color=color, linestyle='dashdot', label='total')

                if D_unc/D > D_ERROR_THRESH:
                    print(f'  total: stopping. D_unc/D = {D_unc/D:.2g}')
                else:

                    print(f'  total:{common.format_val_and_unc(D, D_unc)}, (D_unc/D = {D_unc/D:.2g})')

                    Ds_for_saving.append(D)
                    D_uncs_for_saving.append(D_unc)
                    ks_for_saving.append(k)

            else:
                print(f'  total: stopping. covariance not estimated')

            # new fit
            # new_func = lambda t, a, b, c : a*t**2 + b*t + c
            # # log_func = lambda t, D: np.log10( func(t, D) )
            # # log_unc = lambda x, dx : 0.5 * np.log((x+dx)/(x-dx))
            # # log_f_unc  = log_unc(f, f_unc )
            # if (~f_bad).sum() > 2:
            #     new_popt, new__pcov = scipy.optimize.curve_fit(new_func, t[~f_bad], f[~f_bad])
                
            #     if display:
            #         ax      .plot(t_th, new_func(t_th, *new_popt), color='red', alpha=0.5)
            #         ax_short.plot(t_th, new_func(t_th, *new_popt), color='red', alpha=0.5)

            f_points_short = (~f_bad) & (t < 8)
            f_points_long  = (~f_bad) & (t > 100)

                
            ###### short fit
            if f_points_short.sum() > 2:
                f_popt_short, f_pcov_short = scipy.optimize.curve_fit(
                        log_func, t[f_points_short], np.log10(f[f_points_short]),
                        p0=p0#, sigma=log_f_unc[f_points_short], absolute_sigma=True
                    )
                # print('  short: we had to disable sigma fitting')
                if FIT == EXP_TIMES_CONST_FIT:
                    print(f'  short c = {f_popt_short[1]:.3g} pm {np.sqrt(f_pcov_short[1, 1])}')
                    
                D_short = f_popt_short[0]
                D_unc_short = np.sqrt(f_pcov_short)[0][0]

                if np.isfinite(D_unc_short):
                    
                    if display:
                        ax      .plot(t_th, func(t_th, *f_popt_short), color=color, linestyle='dotted', label=f'short fit {print_D(D_short)}')
                        ax_short.plot(t_th, func(t_th, *f_popt_short), color=color, linestyle='dotted', label=f'short fit {print_D(D_short)}')
                        log_ax  .plot(t_th, func(t_th, *f_popt_short), color=color, linestyle='dotted', label=f'short fit {print_D(D_short)}')
                        D_ax.hlines(D_short,  t.min(), t.max(), color=color, linestyle='dotted', label='short')

                    if D_unc_short/D_short > D_ERROR_THRESH:
                        print(f'  short: stopping. D_unc/D = {D_unc_short/D_short:.2g}')
                    
                    else:
                        print(f'  short: {common.format_val_and_unc(D_short, D_unc_short)}  (D_unc/D = {D_unc_short/D_short:.2g})')
                
                        Ds_for_saving_short    .append(D_short)
                        D_uncs_for_saving_short.append(D_unc_short)
                        ks_for_saving_short    .append(k)
                        
                else:
                    print('  short: stopping. covariance not estimated')

                ###### short linear fit
                if DO_SHORT_LINEAR_FIT:
                    linear_func = lambda t, D : 1 - D*k**2*t
                    f_popt_short_linear, f_pcov_short_linear = scipy.optimize.curve_fit(
                            linear_func, t[f_points_short], f[f_points_short],
                            p0=p0#, sigma=log_f_unc[f_points_short], absolute_sigma=True
                        )
                    D_short_linear = f_popt_short_linear[0]
                    D_unc_short_linear = np.sqrt(f_pcov_short_linear)[0][0]

                    if np.isfinite(D_unc_short_linear):
                        if display:
                            # ax.plot(t_th, linear_func(t_th, *f_popt_short_linear), color=color, linestyle='dotted', label='short')
                            ax_short.plot(t_th, linear_func(t_th, *f_popt_short_linear), color='tab:red', linestyle='dotted', label=f'short lin {print_D(D_short_linear)}')
                            ax      .plot(t_th, linear_func(t_th, *f_popt_short_linear), color='tab:red', linestyle='dotted', label=f'short lin {print_D(D_short_linear)}')
                            log_ax  .plot(t_th, linear_func(t_th, *f_popt_short_linear), color='tab:red', linestyle='dotted', label=f'short lin {print_D(D_short_linear)}')
                            D_ax.hlines(f_popt_short_linear[0],  t.min(), t.max(), color='tab:red', linestyle='dotted', label='short lin')

                        if D_unc_short_linear/D_short_linear > D_ERROR_THRESH:
                            print(f'  short linear: stopping. D_unc/D = {D_unc_short_linear/D_short_linear:.2g}')

                        else:
                            print(f'  short linear: {common.format_val_and_unc(D_short_linear, D_unc_short_linear)} (D_unc/D = {D_unc_short_linear/D_short_linear:.2g})')
                            # Ds_for_saving_short    .append(f_popt_short_linear[0])
                            # D_uncs_for_saving_short.append(np.sqrt(f_pcov_short_linear)[0][0])
                            # ks_for_saving_short    .append(k)
                            
                    else:
                        print('  short linear: stopping. covariance not estimated')
            else:
                print('  short: not attempting')

            ######## long fit
            if f_points_long.sum() > 2:
                f_popt_long, f_pcov_long = scipy.optimize.curve_fit(
                        log_func, t[f_points_long], np.log10(f[f_points_long,]), 
                        # p0=p0, sigma=log_f_unc[f_points_long], absolute_sigma=True
                    )
                if FIT == EXP_TIMES_CONST_FIT:
                    print(f'  long c = {f_popt_long[1]:.3g} pm {np.sqrt(f_pcov_long[1, 1])}')
                
                print('  long:  disabled sigma fitting')
                    
                D_long = f_popt_long[0]
                D_unc_long = np.sqrt(f_pcov_long)[0][0]

                if np.isfinite(D_unc_long):
                    
                    if D_unc_long/D_long > D_ERROR_THRESH:
                        # print(f'  short linear: D_unc/D = {D_unc_long/D_long:.2g}')
                        pass
                    else:
                    
                        if display:
                            ax    .plot(t_th, func(t_th, *f_popt_long), color=color, linestyle='dashed', label=f'long {print_D(D_long)}')
                            log_ax.plot(t_th, func(t_th, *f_popt_long), color=color, linestyle='dashed', label=f'long {print_D(D_long)}')
                        
                        Ds_for_saving_long    .append(D_long)
                        D_uncs_for_saving_long.append(D_unc_long)
                        ks_for_saving_long    .append(k)

                        if display:
                            D_ax.hlines(f_popt_long[0],  t.min(), t.max(), color=color, linestyle='dashed', label='long')
                else:
                    print('  long:  stopping. covariance not estimated')
            else:
                print('  long:  not attempting')


            ######## two stage fit
            if np.sum(~f_bad) < 4:
                print('  twostage: not enough data')
            else:
                func_twostage = lambda t, D1, D2, A1, A2 : A1 * np.exp(-t * k**2 * D1) + A2 * np.exp(-t * k**2 * D2)
                log_func_twostage = lambda t, D1, D2, A1, A2: np.log10( func_twostage(t, D1, D2, A1, A2) )
                # func_twostage = lambda t, D1, D2, A1: A1 * np.exp(-t * k**2 * D1) + (1-A1) * np.exp(-t * k**2 * D2)
                # log_func_twostage = lambda t, D1, D2, A1: np.log10( func_twostage(t, D1, D2, A1) )
                # func_twostage = lambda t, D2, A1: A1 * np.exp(-t * k**2 * D_short) + (1-A1) * np.exp(-t * k**2 * D2)
                # log_func_twostage = lambda t, D2, A1: np.log10( func_twostage(t, D2, A1) )

                try:

                    if D0:
                        p0 = [0.5*D0, 20*D0, 0.9, 0.1]
                    else:
                        p0 = None
                    # p0 = [D_short, 20*D0, 0.1]
                    # p0 = [20*D0, 0.1]
                    f_popt_twostage, f_pcov_twostage = scipy.optimize.curve_fit(
                            log_func_twostage, t[~f_bad], np.log10(f[~f_bad]), p0=p0
                            # p0=p0, sigma=log_f_unc[f_points_long], absolute_sigma=True
                        )
                except RuntimeError as err:
                    print('  twostage:', err)
                
                # D = f_popt[0]
                # D_unc = np.sqrt(f_pcov_twostage)[0][0]
                else:
                    if display and SHOW_TWOSTAGE_FIT:
                        color = 'tab:green'
                        alpha = 1 if np.isfinite(np.sqrt(f_pcov_twostage)[0,0]) else 0.5
                        ax      .plot(t_th, func_twostage(t_th, *f_popt_twostage), color=color, alpha=alpha, linestyle='dashed', label=f'twostage {print_D(f_popt_twostage[0])}, {print_D(f_popt_twostage[1])}')
                        ax_short.plot(t_th, func_twostage(t_th, *f_popt_twostage), color=color, alpha=alpha, linestyle='dashed', label=f'twostage {print_D(f_popt_twostage[0])}, {print_D(f_popt_twostage[1])}')
                        log_ax  .plot(t_th, func_twostage(t_th, *f_popt_twostage), color=color, alpha=alpha, linestyle='dashed', label=f'twostage {print_D(f_popt_twostage[0])}, {print_D(f_popt_twostage[1])}')
                            # ax.plot(t_th, func(t_th, *Fs_popt), color='tab:orange',   linestyle='dotted')
                        # D_ax.hlines(f_popt[0], t.min(), t.max(), color=color, linestyle='dashdot', label='total')

                    # if D_unc/D > D_ERROR_THRESH:
                    #     print(f'  total: stopping. D_unc/D = {D_unc/D:.2g}')
                    # else:
                    if np.isfinite(np.sqrt(f_pcov_twostage)[0,0]):
                        print(f'  twostage: D1={f_popt_twostage[0]:.2g}, D2={f_popt_twostage[1]:.2g}, A1={f_popt_twostage[2]:.2g}, A2={f_popt_twostage[3]:.2g}')
                        # print(f'  twostage: D1={f_popt_twostage[0]:.2g}, D2={f_popt_twostage[1]:.2g}, A1={f_popt_twostage[2]:.2g}')
                        # print(f'  twostage: D2={f_popt_twostage[0]:.2g}, A1={f_popt_twostage[1]:.2g}')

                        # print(f'  total:{common.format_val_and_unc(D, D_unc)}, (D_unc/D = {D_unc/D:.2g})')

                        Ds_for_saving_D1.append(f_popt_twostage[0])
                        D_uncs_for_saving_D1.append(np.sqrt(f_pcov_twostage[0, 0]))
                        ks_for_saving_D1.append(k)
                        
                        Ds_for_saving_D2.append(f_popt_twostage[1])
                        D_uncs_for_saving_D2.append(np.sqrt(f_pcov_twostage[1, 1]))
                        ks_for_saving_D2.append(k)

                    else:
                        print('  twostage: covariance failed')



        ##### displaying

        D      = -1/(k**2 * t ) * np.log(f)
        D_unc  =  1/(k**2 * t ) / np.sqrt(f**2) * f_unc # the sqrt(**2) is needed to prevent negative errors
            
        D2     = -1/k**2 * np.gradient(np.log(f), t)
        D2_unc = np.abs( 1/k**2 * np.gradient(1/f, t) * f_unc )
            
        
        if display:
            D_ax.scatter(t [~f_bad  ], D [~f_bad  ], color=color, label=Ftype, s=6)
            D_ax.scatter(t [ f_bad  ], D [ f_bad  ], color=color,   alpha=0.2, s=6)
            D_ax.scatter(t [~f_bad  ], D2[~f_bad  ], color=colors[type_index+1], label=Ftype, s=6)
            D_ax.scatter(t [ f_bad  ], D2[ f_bad  ], color=colors[type_index+1],   alpha=0.2, s=6)
                # D_ax.errorbar(t , D , yerr=D_unc , fmt='', color=color, alpha=0.2, linestyle='none')
            D_ax.errorbar(t , D , yerr=D_unc , fmt='', color=color, alpha=errorbar_alpha, linestyle='none')
            D_ax.errorbar(t , D2, yerr=D2_unc, fmt='', color=colors[type_index+1], alpha=errorbar_alpha, linestyle='none')
                # D_ax.errorbar(t , D , yerr=D_unc , fmt='', color=color, alpha=0.2, linestyle='none')
                
            D_ax.semilogx()
            log_ax.semilogx()
                # D_ax.hlines(Fs_popt[0], t.min(), t.max(), color='tab:orange', linestyle='dotted')

            # if f_points_short.sum() > 2:
            #     pass
            # else:
            #     print(f'not enough points at k={k:.2f}')

            # if f_points_long.sum() > 2:
            #     pass

                # D_ax.set_ylim(np.nanmin(D), np.nanmax(D))
                
            if file == 'alice0.02':
                D_ax.set_ylim(0, 0.0416*1.6)
                # if file == 'alice0.34':
                #     D_ax.set_ylim(0, 0.031*2)
            if file == 'alice0.66':
                D_ax.set_ylim(0, 0.0175*1.6)
                pass
            if file == 'eleanor0.01' or file == 'eleanor0.34':
                    # D_ax.set_ylim(0, 0.08)
                    # D_ax.set_ylim(0, 1)
                    # D_ax.set_ylim(0, np.nanmax(D))
                pass
            if file == 'eleanor0.34':
                pass
                    # D_ax.set_ylim(0, 0.25)
                    
                    # ax.relim() # tell mpl to ignore errorbars when
                    # ax.autoscale_view() # calcing axis limits
                pass
            D_ax.set_ylim(1e-2, 5e0)
            
            # D_long  = {0.34: 0.023, 0.66: 0.006}
            # D_short = {0.34: 0.033, 0.66: 0.018}
            # D_ax.axhline(y=D_short[phi], color='black', linewidth=1, linestyle=':', label=r'D from MSD')
            # D_ax.axhline(y=D_long [phi], color='black', linewidth=1, linestyle=':')
            
            D_ax    .legend(fontsize=6 if not mult else 5)
            ax_short.legend(fontsize=6)
            log_ax  .legend(fontsize=6)

            ax_short.set_title('$f(k, t)$ lin, short time ' + label, fontsize=8)
            D_ax    .set_title('$D$ ' + label,                       fontsize=8)
            log_ax  .set_title('$f(k, t)$ loglog ' + label,          fontsize=8)
            ax      .set_title('$f(k, t)$ lin ' + label,             fontsize=8)
            
            D_ax.semilogy()

            # D0 = 0.04
            # D_ax.hlines(D0, t.min(), t.max(), color='tab:green', linestyle='dotted')

        if Ftype == 'f' and display and DO_SHORT_LINEAR_FIT:
            # Segre-Pusey rescaling (Segre & Pusey 1996, p772)
            # D2 is D(Q, t) in the paper
            Ds = f_popt_short_linear[0] # Ds(Q) in the paper
            Ds_unc = np.sqrt(f_pcov_short_linear[0, 0])

            low_k = 0.5
            high_k = 10
                
            if low_k <= k*particle_diameter <= high_k or False:
                color = matplotlib.cm.afmhot(np.interp(k*particle_diameter, (low_k, high_k), (0.2, 0.8)))
            
                if SHOW_SEGRE_PUSEY_RESCALING_AXIS:
                    sp_ax = sp_axs[0]
                    sp_q = np.log(f) / (Ds * k**2)
                    
                    sp_ax.plot(t[~f_bad], sp_q[~f_bad], label=k_label, color=color)
                    sp_ax.plot(t[ f_bad], sp_q[ f_bad],                color=color, alpha=0.2)
    
            if SHOW_SEGRE_PUSEY_RESCALING_AXIS:
                Ds_ax = sp_axs[1]

                # Ds_ax.errorbar(k*particle_diameter, Ds, yerr=Ds_unc, marker='o', color='tab:blue')#, label=fr'$k\sigma={k*particle_diameter:.1f}$', color=color)
                Ds_ax.errorbar(k*particle_diameter, D0/Ds, yerr=D0*Ds_unc/Ds**2, marker='o', color='tab:blue')#, label=fr'$k\sigma={k*particle_diameter:.1f}$', color=color)
                    # sp_ax.plot(t[ f_bad], D2[ f_bad]/Ds,                     color=color, alpha=0.2)
    
        if display:
            
            ax.legend(fontsize=6 if not mult else 5)


    if Ftype == 'f' and SHOW_SEGRE_PUSEY_RESCALING_AXIS:
        # sp_ax.semilogx()
        sp_ax.legend()
        sp_ax.set_ylim(-20, 0)
        sp_ax.set_xlim(0, 20)

        sp_ax.set_xlabel('$t$')
        sp_ax.set_ylabel('$\ln(f) / D_S(k) k^2$')
        
        Ds_ax.set_xlabel('$q\sigma$')
        Ds_ax.set_ylabel('$D_0/D_S$')
        # Ds_ax.set_xlim(0, 15)
        # Ds_ax.set_ylim(0, 5)

        common.save_fig(sp_fig, f'scattering_functions/figures_png/segre_pusey_rescaling_{file}.png')     

    if do_fits:
        common.save_data(f'visualisation/data/Ds_from_{Ftype}_{file}',
                    Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, ks=ks_for_saving,
                    particle_diameter=particle_diameter,
                    pixel_size=d.get('pixel_size'),
                    channel=d.get('channel'), NAME=d.get('NAME'))
        common.save_data(f'visualisation/data/Ds_from_{Ftype}_short_{file}',
                    Ds=Ds_for_saving_short, D_uncs=D_uncs_for_saving_short, ks=ks_for_saving_short,
                    particle_diameter=particle_diameter,
                    pixel_size=d.get('pixel_size'),
                    channel=d.get('channel'), NAME=d.get('NAME'))
        common.save_data(f'visualisation/data/Ds_from_{Ftype}_long_{file}',
                    Ds=Ds_for_saving_long, D_uncs=D_uncs_for_saving_long, ks=ks_for_saving_long,
                    particle_diameter=particle_diameter,
                    pixel_size=d.get('pixel_size'),
                    channel=d.get('channel'), NAME=d.get('NAME'))
        common.save_data(f'visualisation/data/Ds_from_{Ftype}_first_{file}',
                    Ds=Ds_for_saving_first, D_uncs=D_uncs_for_saving_first, ks=ks_for_saving_first,
                    particle_diameter=particle_diameter,
                    pixel_size=d.get('pixel_size'),
                    channel=d.get('channel'), NAME=d.get('NAME'))
        common.save_data(f'visualisation/data/Ds_from_{Ftype}_D1_{file}',
                    Ds=Ds_for_saving_D1, D_uncs=D_uncs_for_saving_D1, ks=ks_for_saving_D1,
                    particle_diameter=particle_diameter,
                    pixel_size=d.get('pixel_size'),
                    channel=d.get('channel'), NAME=d.get('NAME'))
        common.save_data(f'visualisation/data/Ds_from_{Ftype}_D2_{file}',
                    Ds=Ds_for_saving_D2, D_uncs=D_uncs_for_saving_D2, ks=ks_for_saving_D2,
                    particle_diameter=particle_diameter,
                    pixel_size=d.get('pixel_size'),
                    channel=d.get('channel'), NAME=d.get('NAME'))

if __name__ == '__main__':
    for file in sys.argv[1:]:

        
        num_displayed_ks = 20
        fig, axes = plt.subplots(4, num_displayed_ks-5, figsize=(num_displayed_ks*3, 4*2.8))
        #                         big big big hack: ^^

        for type_index, Ftype in enumerate(['f', 'Fs']):
        # for type_index, Ftype in enumerate(['Fs', 'f', 'DDM']):

            fig = show_single_F_type(file, type_index, Ftype, fig, axes, num_displayed_ks)

                
        plt.suptitle(fr'F or F_s (k, t), {file}')

        common.save_fig(fig, f'scattering_functions/figures_png/Fs_decay_t_{file}.png', dpi=300)
        