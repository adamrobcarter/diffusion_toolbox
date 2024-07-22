import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings, math
import scipy.optimize, scipy.signal
import matplotlib.cm

subplot_i = 0

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)
# target_ks = list(np.logspace(np.log10(0.02), np.log10(8), 20))
# target_ks = list(np.logspace(np.log10(0.02), np.log10(0.1), 25))

def show_single_F_type(file, type_index, Ftype, fig, axes, num_displayed_ks, mult=False):
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

    if Ftype == 'f':
        sp_fig, sp_axs = plt.subplots(2, 1, figsize=(4, 6))

    load = 'F' if Ftype=='f' else Ftype
    d = common.load(f"scattering_functions/data/{load}_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]
    particle_diameter     = d['particle_diameter']

    every_nth_k = int(math.ceil(k_all.shape[1] / num_displayed_ks))
    every_nth_k = max(every_nth_k, 1)

    graph_i = 0

    for k_index in range(k_all.shape[1]):
        # target_k = target_ks[graph_i]
            # k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

        if Ftype == 'DDM':
            ks = k_all
        else:
            ks = k_all[0, :]

        # k_index = np.argmax(ks >= target_k)
        k = ks[k_index]


        # print(f'k: target {target_k:.3f}, real {k:.3f}, index {k_index}, 2pi/k={2*np.pi/k:.1f}um')
        print(f'k {k:.3f}, index {k_index}, 2pi/k={2*np.pi/k:.1f}um')

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

        display = False
        if k_index % every_nth_k == 0:
            display = True

            ax = lin_axes[graph_i]
            ax_short = lin_short_axes[graph_i]
            D_ax = D_axes[graph_i]
            log_ax = log_axes[graph_i]
            log_ax.semilogy()
            # extra_ax = extra_axes[graph_i]

            graph_i += 1

        label = fr"$k={k:.3f}\mathrm{{\mu m}}$ ($L\approx{2*np.pi/k:.1f}\mathrm{{\mu m}}$, $k\sigma={k*particle_diameter:.2f}$)"

        if common.nanfrac(f) == 1:
            print(f'all nan at k={k:.1f} (i={graph_i})')
            continue


        if Ftype == 'f':
            noise_thresh = 1e-2 # for eleanorlong!!
            time_thresh = 200
        elif Ftype == 'Fs':
            noise_thresh = 1.7e-2
            time_thresh = 400
        elif Ftype == 'DDM':
            noise_thresh = 1e-3
            time_thresh = 400

        f_noise   = f < noise_thresh
        f_toolong = t > time_thresh
            # f_bad   = f_noise | f_toolong # f_toolong should depend on k!!!!
        f_bad   = f_noise
        f_bad[0] = True

        # new noise identification idea
        # first noise point is first point where the gradient is no longer getting more negative
        grad = np.gradient(np.log10(f), np.log10(t))
        peaks, props = scipy.signal.find_peaks(-grad, prominence=0.01)
        # prominance filters out peaks that are just noise. A smoothing filter would probably be better though
        f_bad = np.full(f.shape, False)
        print(np.log10(t[np.array(peaks)]))
        print(props)
        if len(peaks):
            f_bad[peaks[0]:] = True
        else:
            print('  no peaks in -grad found?!')
        # print(f_bad, f_bad.dtype)
        # if display:
        #     extra_ax.scatter(np.log10(t), grad, s=3)
        # extra_axes[graph_i].semilogx()

        f_bad[f<0] = True

        if display:
            ax.plot(t[~f_bad], f[~f_bad], color=colors[type_index], linestyle='', label=f'{Ftype} {file}', marker='.')
            ax.errorbar(t[~f_bad], f[~f_bad], yerr=f_unc[~f_bad], color=colors[type_index], linestyle='', alpha=0.2, marker='none')
            ax.errorbar(t[ f_bad], f[ f_bad], yerr=f_unc[ f_bad], color=colors[type_index], linestyle='', alpha=0.2,        marker='.')

        
            log_ax.plot(t[~f_bad], f[~f_bad], color=colors[type_index], linestyle='', label=f'{Ftype} {file}', marker='.')
            log_ax.errorbar(t[~f_bad], f[~f_bad], yerr=f_unc[~f_bad], color=colors[type_index], linestyle='', alpha=0.2, marker='none')
            log_ax.errorbar(t[ f_bad], f[ f_bad], yerr=f_unc[ f_bad], color=colors[type_index], linestyle='', alpha=0.2,        marker='.')
            log_ax.set_ylim(f.min()/1.2, f.max()*1.2)

            # ax.errorbar(t[~f_bad], f[~f_bad], yerr=0,           color=colors[type_index], linestyle='', label=f'{Ftype} {file}', marker='.')
            # ax.errorbar(t[ f_bad], f[ f_bad], yerr=0,           color=colors[type_index], linestyle='', alpha=0.2,        marker='.')
            # ax.scatter(t[~f_bad], f[~f_bad], color=colors[type_index], label=f'{Ftype} {file}', s=6)
            # ax.scatter(t[ f_bad], f[ f_bad], color=colors[type_index], alpha=0.3,        s=6)
            
            # print(f'{Ftype} k={k:.2f} avg err = ')

            ax.legend(fontsize=7 if not mult else 5)
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
            end_plot_time = 1/(k+0.1) * 500
            end_plot_y = f[np.argmax(t > end_plot_time)] * 0.8
            if end_plot_y > 1: end_plot_y = 0.9
            if end_plot_y < 1e-4: end_plot_y = 1e-4
            ax.set_ylim(end_plot_y , 1.01)
                # ax.set_xlim(0, 1/k**2 * 100)
            ax.set_xlim(0, end_plot_time)

            # ax.set_ylim(min(max(func(t[np.argmax(f_bad[1:])], *f_popt), 1e-3), 9.9e-1), 1.01)
            # ax.set_ylim(1e-3, 1.1)



        negative = f <= 0
        if negative.sum() > 0:
            print(f'  negative: {negative.sum()/negative.size:.2f}')

        # fits
        # we fit to log(f(k, t)) because otherwise points near one make a much
        # larger impact to the fit than points near zero
        func = lambda t, D : np.exp(-t * k**2 * D)
        log_func = lambda t, D: np.log10( func(t, D) )
        log_unc = lambda x, dx : 0.5 * np.log((x+dx)/(x-dx))
        log_f_unc  = log_unc(f, f_unc )
            
        p0 = [0.05]
        p0 = None

        if (~f_bad).sum() == 0:
            print('  no good data, skipping rest')
            continue
            
        f_popt, f_pcov = scipy.optimize.curve_fit(
                log_func, t[~f_bad], np.log10(f[~f_bad]),
                # p0=p0, sigma=log_f_unc[~f_bad], absolute_sigma=True
            )
        print('  total: disabled sigma fitting')
            # Fs_popt, Fs_pcov = scipy.optimize.curve_fit(log_func, t2[~Fs_bad], np.log10(Fs[~Fs_bad]), sigma=log_Fs_unc[~Fs_bad], absolute_sigma=True)
        t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))

        if np.isfinite(np.sqrt(f_pcov)[0,0]):
            print('  total:', common.format_val_and_unc(f_popt[0], np.sqrt(f_pcov)[0][0]))
            
            if display:
                ax.plot(t_th, func(t_th, *f_popt), color=colors[type_index], linestyle='dashdot', label='total')
                    # ax.plot(t_th, func(t_th, *Fs_popt), color='tab:orange',   linestyle='dotted')

            Ds_for_saving.append(f_popt[0])
            D_uncs_for_saving.append(np.sqrt(f_pcov)[0][0])
            ks_for_saving.append(k)

            if display:
                D_ax.hlines(f_popt[0], t.min(), t.max(), color=colors[type_index], linestyle='dashdot', label='total')

        else:
            print(f'  total: covariance not estimated')

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


        if display:
            SHORT_PLOT_MAX_X = 20
            ax_short.errorbar(t, f, yerr=f_unc, color=colors[type_index], linestyle='', alpha=0.2)
            ax_short.errorbar(t, f, yerr=0, color=colors[type_index], linestyle='', marker='.')
            # ax_short.errorbar(t, f, yerr=0, color=colors[type_index], linestyle='', marker='.')
            # ax_short.semilogy()
            ax_short.set_xlim(0, SHORT_PLOT_MAX_X)
            ax_short.set_ylim(f[0:SHORT_PLOT_MAX_X+1].min(), f[0:SHORT_PLOT_MAX_X+1].max())
            
        if f_points_short.sum() > 2:
            f_popt_short, f_pcov_short = scipy.optimize.curve_fit(
                    log_func, t[f_points_short], np.log10(f[f_points_short]),
                    p0=p0, sigma=log_f_unc[f_points_short], absolute_sigma=True
                )
            # print('  short: we had to disable sigma fitting')
                
            if np.isfinite(np.sqrt(f_pcov_short)[0][0]):
                
                if display:
                    ax.plot(t_th, func(t_th, *f_popt_short), color=colors[type_index], linestyle='dotted', label='short')
                    ax_short.plot(t_th, func(t_th, *f_popt_short), color=colors[type_index], linestyle='dotted', label='short')
        
                print(f'  short: {common.format_val_and_unc(f_popt_short[0], np.sqrt(f_pcov_short)[0][0])}')
                Ds_for_saving_short    .append(f_popt_short[0])
                D_uncs_for_saving_short.append(np.sqrt(f_pcov_short)[0][0])
                ks_for_saving_short    .append(k)
                    
                if display:
                    D_ax.hlines(f_popt_short[0],  t.min(), t.max(), color=colors[type_index], linestyle='dotted', label='short')
            else:
                print('  short: covariance not estimated')

                
            linear_func = lambda t, D : 1 - D*k**2*t
            f_popt_short_linear, f_pcov_short_linear = scipy.optimize.curve_fit(
                    linear_func, t[f_points_short], f[f_points_short],
                    p0=p0#, sigma=log_f_unc[f_points_short], absolute_sigma=True
                )
                
            if np.isfinite(np.sqrt(f_pcov_short_linear)[0][0]):
                
                if display:
                    # ax.plot(t_th, linear_func(t_th, *f_popt_short_linear), color=colors[type_index], linestyle='dotted', label='short')
                    ax_short.plot(t_th, linear_func(t_th, *f_popt_short_linear), color='tab:red', linestyle='dotted', label='short lin')
        
                print(f'  short linear: {common.format_val_and_unc(f_popt_short_linear[0], np.sqrt(f_pcov_short_linear)[0][0])}')
                # Ds_for_saving_short    .append(f_popt_short_linear[0])
                # D_uncs_for_saving_short.append(np.sqrt(f_pcov_short_linear)[0][0])
                # ks_for_saving_short    .append(k)
                    
                if display:
                    D_ax.hlines(f_popt_short_linear[0],  t.min(), t.max(), color='tab:red', linestyle='dotted', label='short lin')
            else:
                print('  short linear: covariance not estimated')
        else:
            print('  short: not attempting')

        if f_points_long.sum() > 2:
            f_popt_long, f_pcov_long = scipy.optimize.curve_fit(
                    log_func, t[f_points_long], np.log10(f[f_points_long,]), 
                    # p0=p0, sigma=log_f_unc[f_points_long], absolute_sigma=True
                )
            print('  long:  disabled sigma fitting')
                
            if np.isfinite(np.sqrt(f_pcov_long)[0][0]):
                
                if display:
                    ax.plot(t_th, func(t_th, *f_popt_long), color=colors[type_index], linestyle='dashed', label='long')
                
                Ds_for_saving_long    .append(f_popt_long[0])
                D_uncs_for_saving_long.append(np.sqrt(f_pcov_long)[0][0])
                ks_for_saving_long    .append(k)

                if display:
                    D_ax.hlines(f_popt_long[0],  t.min(), t.max(), color=colors[type_index], linestyle='dashed', label='long')
            else:
                print('  long:  covariance not estimated')
        else:
            print('  long:  not attempting')


            
        D      = -1/(k**2 * t ) * np.log(f)
        D_unc  =  1/(k**2 * t ) / np.sqrt(f**2) * f_unc # the sqrt(**2) is needed to prevent negative errors
            
        D2     = -1/k**2 * np.gradient(np.log(f), t)
        D2_unc = np.abs( 1/k**2 * np.gradient(1/f, t) * f_unc )
            
        
        if display:
            D_ax.scatter(t [~f_bad  ], D [~f_bad  ], color=colors[type_index], label=Ftype, s=6)
            D_ax.scatter(t [ f_bad  ], D [ f_bad  ], color=colors[type_index],   alpha=0.2, s=6)
            D_ax.scatter(t [~f_bad  ], D2[~f_bad  ], color=colors[type_index+1], label=Ftype, s=6)
            D_ax.scatter(t [ f_bad  ], D2[ f_bad  ], color=colors[type_index+1],   alpha=0.2, s=6)
                # D_ax.errorbar(t , D , yerr=D_unc , fmt='', color=colors[type_index], alpha=0.2, linestyle='none')
            D_ax.errorbar(t , D , yerr=D_unc , fmt='', color=colors[type_index], alpha=0.2, linestyle='none')
            D_ax.errorbar(t , D2, yerr=D2_unc, fmt='', color=colors[type_index+1], alpha=0.2, linestyle='none')
                # D_ax.errorbar(t , D , yerr=D_unc , fmt='', color=colors[type_index], alpha=0.2, linestyle='none')
                
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
            D_ax.set_ylim(1e-2, 1e0)
            
            # D_long  = {0.34: 0.023, 0.66: 0.006}
            # D_short = {0.34: 0.033, 0.66: 0.018}
            # D_ax.axhline(y=D_short[phi], color='black', linewidth=1, linestyle=':', label=r'D from MSD')
            # D_ax.axhline(y=D_long [phi], color='black', linewidth=1, linestyle=':')
            
            D_ax.legend(fontsize=7 if not mult else 5)
            D_ax.semilogy()

            D0 = 0.04
            D_ax.hlines(D0, t.min(), t.max(), color='tab:green', linestyle='dotted')

        if Ftype == 'f' and display:
            # Segre-Pusey rescaling (Segre & Pusey 1996, p772)
            # D2 is D(Q, t) in the paper
            Ds = f_popt_short_linear[0] # Ds(Q) in the paper
            Ds_unc = np.sqrt(f_pcov_short_linear[0, 0])

            low_k = 0.5
            high_k = 10
                
            if low_k <= k*particle_diameter <= high_k or False:
                color = matplotlib.cm.afmhot(np.interp(k*particle_diameter, (low_k, high_k), (0.2, 0.8)))
            
                sp_ax = sp_axs[0]
                sp_q = np.log(f) / (Ds * k**2)
                
                sp_ax.plot(t[~f_bad], sp_q[~f_bad], label=fr'$k\sigma={k*particle_diameter:.1f}$', color=color)
                sp_ax.plot(t[ f_bad], sp_q[ f_bad],                     color=color, alpha=0.2)
    
            Ds_ax = sp_axs[1]

            # Ds_ax.errorbar(k*particle_diameter, Ds, yerr=Ds_unc, marker='o', color='tab:blue')#, label=fr'$k\sigma={k*particle_diameter:.1f}$', color=color)
            Ds_ax.errorbar(k*particle_diameter, D0/Ds, yerr=D0*Ds_unc/Ds**2, marker='o', color='tab:blue')#, label=fr'$k\sigma={k*particle_diameter:.1f}$', color=color)
                # sp_ax.plot(t[ f_bad], D2[ f_bad]/Ds,                     color=color, alpha=0.2)
    

    if Ftype == 'f':
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

    common.save_data(f'visualisation/data/Ds_from_{Ftype}_{file}',
                Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, ks=ks_for_saving,
                particle_diameter=particle_diameter,
                pixel_size=d.get('pixel_size'))
    common.save_data(f'visualisation/data/Ds_from_{Ftype}_short_{file}',
                Ds=Ds_for_saving_short, D_uncs=D_uncs_for_saving_short, ks=ks_for_saving_short,
                particle_diameter=particle_diameter,
                pixel_size=d.get('pixel_size'))
    common.save_data(f'visualisation/data/Ds_from_{Ftype}_long_{file}',
                Ds=Ds_for_saving_long, D_uncs=D_uncs_for_saving_long, ks=ks_for_saving_long,
                particle_diameter=particle_diameter,
                pixel_size=d.get('pixel_size'))

if __name__ == '__main__':
    for file in sys.argv[1:]:

        
        num_displayed_ks = 20
        fig, axes = plt.subplots(4, num_displayed_ks, figsize=(num_displayed_ks*3, 4*2.8))
        
        for type_index, Ftype in enumerate(['f', 'Fs']):
        # for type_index, Ftype in enumerate(['Fs', 'f', 'DDM']):

            fig = show_single_F_type(file, type_index, Ftype, fig, axes, num_displayed_ks)

                
        plt.suptitle(fr'F or F_s (k, t), {file}')

        common.save_fig(fig, f'scattering_functions/figures_png/Fs_decay_t_{file}.png', dpi=300)
        