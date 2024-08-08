import countoscope_theory.nmsd
import countoscope_theory.structure_factor
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
# import sDFT_interactions
import matplotlib.cm

import countoscope_theory

PRESENT_SMALL = False
LABELS_ON_PLOT = False
SHOW_THEORY_FIT = False
SHOW_PLATEAUS_THEORY = False
SHOW_VARIANCE = False
SHOW_MEAN = False
SHOW_PLATEAUS_OBS = False
SHOW_PLATEAU_OBS_AREA = False
SHOW_SHORT_TIME_FIT = False

figsize = (6, 4.5)
if PRESENT_SMALL:
    figsize = (4.5, 4)
    figsize = (3.5, 3.2)

collapse_x = True
collapse_y = True
# collapse_x = False
# collapse_y = False

def get_plateau(nmsd, file):
                
    if file == 'eleanorlong': # hacky, pls don't do this
        start_index = -70000
        end_index   = -20000
    else: # used to be -300, -100
        start_index = -600
        end_index   = -400

    used_data = nmsd[start_index:end_index] # used to be -300:-100, we could do with a more inteligent method (use the gradient (smoothed?))
    # used_data = nmsd[-300:-100]
    return used_data.mean(), used_data.std()

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'counted_'):

        D0_from_fits     = [{}, {}]
        D0_unc_from_fits = [{}, {}]

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
        N_mean       = data['N_mean']
        N_var        = data['N_var']
        N_var_mod    = data['N_var_mod']
        num_of_boxes = data['num_boxes']
        sep_sizes    = data['sep_sizes']

        num_timesteps = N2_mean.shape[1]
        num_boxes     = N2_mean.shape[0]
        t_all = np.arange(0, num_timesteps) * time_step

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


        # for box_size_index, L in enumerate(box_sizes):
        for box_size_index, L in list(enumerate(box_sizes)):
        # for L in [2**e for e in range(-2, 7)]:

            L   = box_sizes[box_size_index]
            sep = sep_sizes[box_size_index]

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

            #, r^2={r2:.2f}
            label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
            # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

            # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

            # color = matplotlib.cm.afmhot((box_size_index+2)/(len(box_sizes)+7))
            color =  matplotlib.cm.afmhot(np.interp(box_size_index, (0, len(box_sizes)), (0.2, 0.75)))
                    
            if SHOW_MEAN:
                ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color=color, linewidth=1, linestyle='dashdot', label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
            if SHOW_VARIANCE:
                ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$2\mathrm{Var}(N)$' if box_size_index==0 else None)
                ax.hlines(2*N_var_mod[box_size_index], t.min(), t.max(), linestyles='dotted', color='grey', linewidth=1, label=r'$2\mathrm{Var}(N)$' if box_size_index==0 else None)


            # N2_theory = 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2)**2) # countoscope eq. 2
            # N2_theory_lowtime = 4 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t_theory) * (L + L_2) / (L * L_2)
            # ax.plot(t_theory[:LOWTIME_FIT_END], N2_theory_lowtime[:LOWTIME_FIT_END], linestyle='dashed', linewidth=1, color='black', label='sFDT (no inter.) low time' if box_size_index==0 else None)

            # p1, p2 = plateaus.calc_plateaus_for_L(sigma, phi, L)
            # ax.hlines(p1, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1, label='plateaus')
            # ax.hlines(p2, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1)

            # computed theory interactions
            # D0 = { # countoscope paper, table 1
            #     'alice0.02': 0.0416,
            #     'alice0.02_overlapped': 0.0416,
            #     'alice0.34': 0.0310,
            #     'alice0.66': 0.0175
            # }[file]
            # N2_theory_interactions = 2 * N_mean[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
            # ax.plot(t_theory, N2_theory_interactions, color='black', linewidth=1, label='sFDT (w/ inter.)' if box_size_index==0 else None)

            # fit to whole thing
            if depth_of_field:
                N2_theory = lambda t, D, N: common.N2_nointer_3D(t, D, N, L, L, depth_of_field)
                type_of_fit = 'sDFT fit (no inter, 3D)'
            else:
                N2_theory = lambda t, D : countoscope_theory.nmsd.inter_2d(t, D, N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
                type_of_fit = 'sDFT fit (w/ inter.)'
            log_N2_theory = lambda t, *args : np.log(N2_theory(t, *args)) # we fit to log otherwise the smaller points make less impact to the fit
            
            fitting_points = common.exponential_integers(1, t.size//2)
            # p0 = (0.05, N_mean[box_size_index])
            p0 = [0.05]
            popt, pcov = scipy.optimize.curve_fit(log_N2_theory, t[fitting_points], np.log(delta_N_sq[fitting_points]), p0=p0, maxfev=2000)
            # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
            D_from_fit = popt[0]
            D_from_fit_unc = np.sqrt(pcov[0][0])
            
            N2_theory_points = N2_theory(t_theory, *popt)
            
            if collapse_y:
                rescale = N_mean[box_size_index]
                # rescale = N_var[box_size_index]
                delta_N_sq       /= rescale
                delta_N_sq_err   /= rescale
                N2_theory_points /= rescale
            if collapse_x:
                t /= L**2
                t_theory /= L**2
            
            
            Ds_for_saving.append(D_from_fit)
            D_uncs_for_saving.append(D_from_fit_unc)
            Ls_for_saving.append(L)
            
            info = fr'L = {L/sigma:.1f}σ, D_fit = {common.format_val_and_unc(D_from_fit, D_from_fit_unc, 2)} um^2/s'
            # ±{np.sqrt(pcov[0][0]):.3f}$'
            # print(info)

            
            # linear fit to start
            fit_end = 6
            fit_func_2 = lambda t, D : 8 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(t * D / L**2)
            popt, pcov = scipy.optimize.curve_fit(fit_func_2, t[1:fit_end], delta_N_sq[1:fit_end])
            D_from_shorttime = popt[0]
            D_unc_from_shorttime = np.sqrt(pcov)[0, 0]
            # print(f'L={L}um, D_unc/D={D_unc_from_shorttime/D_from_shorttime:.3f}')
            if D_unc_from_shorttime/D_from_shorttime > 0.03:
                pass
            else:
                Ds_shorttime_for_saving.append(D_from_shorttime)
                D_uncs_shorttime_for_saving.append(D_unc_from_shorttime)
                Ls_shorttime_for_saving.append(L)
            

            if len(box_sizes) <= 10:
                display = True
            else:
                display = box_size_index % (len(box_sizes) // 10) == 0

            # display = box_size_index == 20

            # display = True
            if display:

                if not collapse_y and SHOW_THEORY_FIT:
                    ax.plot(t_theory[1:], N2_theory_points[1:], color='black', linewidth=1, label=type_of_fit if box_size_index==0 else None)
                # label += fr', $D_\mathrm{{fit}}={popt[0]:.3g}\pm {np.sqrt(pcov[0][0]):.3g} \mathrm{{\mu m^2/s}}$'

                if collapse_x or collapse_y:
                    markersize = 2
                else:
                    if PRESENT_SMALL:
                        markersize = 5
                    else:
                        markersize = 3

                # actual data
                if LABELS_ON_PLOT:
                    label = label='observations' if box_size_index==0 else None
                else:
                    label = f'L={L:.2f}'
                    label += f', s={sep:.2f}'
                label += fr', $D_\mathrm{{short\:fit}}={common.format_val_and_unc(D_from_fit, D_from_fit_unc, 2)} \mathrm{{\mu m^2/s}}$'
                # ±{np.sqrt(pcov[0][0]):.3f}$'
                print(delta_N_sq.size, common.nanfrac(delta_N_sq))
                exp_plot = ax.plot(t[1:], delta_N_sq[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1, color=color)
                # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:]/np.sqrt(num_of_boxes[box_size_index]), label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
                # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            
                if LABELS_ON_PLOT:
                    t_index_for_text = int(t_theory.size // 1.6)
                    angle = np.tan(np.gradient(N2_theory_points, t_theory)[t_index_for_text]) * 180/np.pi
                    # plt.scatter(t_theory[t_index_for_text], N2_theory_points[t_index_for_text])
                    L_label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
                    L_label = rf'$L={L/sigma:.1f}\sigma$'
                    ax.text(t_theory[t_index_for_text+6], N2_theory_points[t_index_for_text+6]*1.4, L_label,
                            horizontalalignment='center', color=color, fontsize=9,
                            transform_rotates_text=True, rotation=angle, rotation_mode='anchor')
                
                # linear fit to start
                if D_unc_from_shorttime/D_from_shorttime > 0.03:
                    print(f'skipping short time fit at L={L}um, D_unc/D={D_unc_from_shorttime/D_from_shorttime:.2f}')
                else:
                    D_ratio = D_from_shorttime/D_from_fit
                    # print(f'D_short / D_fit = {D_ratio:.2f}')
                    if D_ratio > 1.5 or 1/D_ratio > 1.5:
                        print(f'problem! D fit = {common.format_val_and_unc(D_from_fit, D_from_fit_unc, 2)} D shorttime = {common.format_val_and_unc(D_from_shorttime, D_unc_from_shorttime, 2)}')
                    if SHOW_SHORT_TIME_FIT:
                        ax.plot(t[1:fit_end], fit_func_2(t[1:fit_end], *popt), linestyle=':', color='gray')

            if SHOW_PLATEAUS_THEORY:
                ax.hlines(
                    countoscope_theory.nmsd.plateau_inter_2d(N_mean[box_size_index], L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma)),
                    t[0], t[-1], linestyle='dotted', color=color, label='sDFT plateau' if box_size_index==0 else None)
            if SHOW_PLATEAUS_OBS:

                plat, plat_std = get_plateau(N2_mean[box_size_index, :], PLATEAU_OBS_START, PLATEAU_OBS_END)
                ax.hlines(plat, t[0], t[-1], linestyle='dotted', color='white', label='obs plat' if box_size_index==0 else None)
        

        if not collapse_x:
            ax.legend(fontsize=5 if not PRESENT_SMALL else 5, loc='lower left')
        ax.semilogy()
        ax.semilogx()
        xlabel = '$t/L^2$' if collapse_x else '$t$ ($\mathrm{s}$)'
        ylabel = r'$\Delta N^2(t)/ Var(N)$' if collapse_y else r'$\langle \Delta N^2(t) \rangle$ ($\mathrm{\mu m^2}$)'
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

        # common.save_fig(fig, f'/home/acarter/presentations/cmd31/figures/nmsd_one_{file}.pdf', hide_metadata=True)
        common.save_fig(fig, f'box_counting/figures_png/nmsd_{file}.png', dpi=200)

        common.save_data(f'visualisation/data/Ds_from_boxcounting_{file}',
                Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, Ls=Ls_for_saving,
                particle_diameter=sigma)
        common.save_data(f'visualisation/data/Ds_from_boxcounting_shorttime_{file}',
                Ds=Ds_shorttime_for_saving, D_uncs=D_uncs_shorttime_for_saving, Ls=Ls_shorttime_for_saving,
                particle_diameter=sigma)
        