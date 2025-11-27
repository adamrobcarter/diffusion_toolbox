import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
# import sDFT_interactions
import matplotlib.cm

PRESENT_SMALL = False
LABELS_ON_PLOT = False
SHOW_THEORY_FIT = True
SHOW_PLATEAUS = False
SHOW_VARIANCE = False
SHOW_MEAN = False
LOGLOG = False
ZOOM_TO_START = True
# ZOOM_TO_START = False

figsize = (6, 4.5)
if PRESENT_SMALL:
    figsize = (4.5, 4)
    figsize = (3.5, 3.2)

collapse_x = True
collapse_y = True
collapse_x = False
# collapse_y = False

I_0 = 1100 # from Sophie in slack

for file in common.files_from_argv('intensity_counting/data/', 'counted_'):

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    # rescaled_fig, rescaled_axs = plt.subplots(2, 1, figsize=(5, 8), squeeze=False)

    data = common.load(f'intensity_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean        = data['counted_intensity_diffs']
    # N2_std         = data['N2_std']
    time_step      = data['time_step']
    depth_of_field = data.get('depth_of_field')

    box_sizes    = data['box_sizes']
    I_mean       = data['avg_intensities']
    I_var        = data['variances']

    print('I_mean', I_mean)
    # sep_sizes    = data['sep_sizes']

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
        # sep = sep_sizes[box_size_index]

        delta_N_sq     = N2_mean[box_size_index, :]
        # delta_N_sq_err = N2_std [box_size_index, :]
        t = np.copy(t_all)
        t_theory = np.logspace(np.log10(t_all[1] / 2), np.log10(t_all.max()*1), 100)


        
        if collapse_y:
            delta_N_sq       /= I_mean[box_size_index] * I_0
            print(delta_N_sq)
            # delta_N_sq_err   /= N_var[box_size_index]
        if collapse_x:
            t        = t        / L**2
            t_theory = t_theory / L**2

        anomalous = delta_N_sq < 1e-14
        anomalous[0] = False # don't want to remove point t=0 as it could legit be zero
        if np.any(anomalous):
            print(f'found {anomalous.sum()/delta_N_sq.size*100:.3f}% anomalous')
            delta_N_sq     = delta_N_sq    [~anomalous]
            # delta_N_sq_err = delta_N_sq_err[~anomalous]
            t              = t             [~anomalous]
        
        #, r^2={r2:.2f}
        label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
        # label += f', $D={D0:.3f}±{np.sqrt(pcov[0][0]):.3f}$'

        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        # color = matplotlib.cm.afmhot((box_size_index+2)/(len(box_sizes)+7))
        color =  matplotlib.cm.afmhot(np.interp(box_size_index, (0, len(box_sizes)), (0.2, 0.75)))
                
        if SHOW_MEAN:
            ax.hlines(2*I_mean[box_size_index], t.min(), t.max(), color=color, linewidth=1, linestyle='dashdot', label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        if SHOW_VARIANCE:
            ax.hlines(2*I_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$2\mathrm{Var}(N)$' if box_size_index==0 else None)

        
        if (not PRESENT_SMALL and len(box_sizes)<15) or box_size_index % (len(box_sizes) // 10) == 0:
            # only plot sometimes
            
            
            # if collapse_x or collapse_y:
            #     markersize = 2
            # else:
            #     if PRESENT_SMALL:
            #         markersize = 5
            #     else:
            #         markersize = 3
            markersize = 5

            # actual data
            if LABELS_ON_PLOT:
                label = label='observations' if box_size_index==0 else None
            else:
                label = f'L={L:.2f}'
                # label += f', s={sep:.2f}'
            
            print('size, nanfrac', delta_N_sq.size, common.nanfrac(delta_N_sq))
            exp_plot = ax.plot(t[:], delta_N_sq[:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1, color=color)
            # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:]/np.sqrt(num_of_boxes[box_size_index]), label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
            # exp_plot = ax.errorbar(t[1:], delta_N_sq[1:], yerr=delta_N_sq_err[1:], label=label, linestyle='none', marker='o', markersize=markersize, zorder=-1)
        
            # if LABELS_ON_PLOT:
            #     t_index_for_text = int(t_theory.size // 1.6)
            #     angle = np.tan(np.gradient(N2_theory_points, t_theory)[t_index_for_text]) * 180/np.pi
            #     # plt.scatter(t_theory[t_index_for_text], N2_theory_points[t_index_for_text])
            #     L_label = rf'$L={L:.1f}\mathrm{{\mu m}}$'
            #     L_label = rf'$L={L/sigma:.1f}\sigma$'
            #     ax.text(t_theory[t_index_for_text+6], N2_theory_points[t_index_for_text+6]*1.4, L_label,
            #             horizontalalignment='center', color=color, fontsize=9,
            #             transform_rotates_text=True, rotation=angle, rotation_mode='anchor')
            
    # if collapse_y:
    #     ax.plot([1e-3, 1e-2], [1e-1, 1e0], color='blue', label='t^1')
    #     ax.plot([1e-3, 1e-2], [1e-1, np.sqrt(10)*1e-1], color='green', label='t^1/2')

    ax.legend(fontsize=6 if not PRESENT_SMALL else 7, loc='lower right')
    # ax.semilogy()
    # ax.set_ylim(0, 0.02)
    # ax.semilogx()
    xlabel = '$t/L^2$' if collapse_x else r'$t$ ($\mathrm{s}$)'
    ylabel = r'$\Delta I^2(t)/\langle I \rangle$' if collapse_y else r'$\langle \Delta I^2(t) \rangle$ ($\mathrm{\mu m^2}$)'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.grid()
    if LOGLOG:
        ax.semilogx()
        ax.semilogy()

    if ZOOM_TO_START:
        if collapse_x:
            pass
            ax.set_xlim(0, 0.01)
        else:
            ax.set_xlim(0, 2)

        if collapse_y:
            pass
            ax.set_ylim(0, 0.05)
        else:
            ax.set_ylim(0, 100)
    
    title = file
    # title = f'Simulated colloids in RCP spheres\n$\phi={phi:.3f}$'
    # if not np.isnan(phi):
    #     title += f', $\phi_\mathrm{{calc}}={phi:.3f}$'
    # if not np.isnan(sigma):
    #     title += f', $\sigma={sigma:.2f}\mathrm{{\mu m}}$'
    # if sigma_calced := data.get('particle_diameter_calced'):
    #     title += f', $\sigma_\mathrm{{calc}}={sigma_calced:.3f}\mathrm{{\mu m}}$'
        # print('sigma calced hidden from legend')
    if not PRESENT_SMALL:
        ax.set_title(title)

    # common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/boxcounting_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'intensity_counting/figures_png/msd_{file}.png', dpi=200)

    # common.save_data(f'visualisation/data/Ds_from_boxcounting_{file}',
    #          Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, Ls=Ls_for_saving,
    #          particle_diameter=sigma)
    # common.save_data(f'visualisation/data/Ds_from_boxcounting_shorttime_{file}',
    #          Ds=Ds_shorttime_for_saving, D_uncs=D_uncs_shorttime_for_saving, Ls=Ls_shorttime_for_saving,
    #          particle_diameter=sigma)
    