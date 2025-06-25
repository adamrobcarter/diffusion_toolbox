import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.special
import math

FIT_USE_FLOW = True
D_ERROR_THRESH = 0.5 # D_unc/D must be smaller than this to save
SEMILOGX = True
# THEORY_PLOT_DS = [1, 0.1]
THEORY_PLOT_DS = []
DO_FIT = True

def plot_single_graph(ax, t, k, F_D_sq, F_D_sq_unc, k_index, plot_label=None, color=None, do_fit=False, particle_diameter=None,
                      particle_material=None, fit_use_flow=False):
    to_plot = np.full_like(F_D_sq[:, k_index], True, dtype='bool')

    if SEMILOGX:
        to_plot[0] = False
    anomalous = F_D_sq_unc[:, k_index] > F_D_sq[:, k_index]
    # to_plot[anomalous] = False


    if num_removed := (to_plot==False).sum()-1:
        print(f'  removed {num_removed}')


    ax.errorbar(t[to_plot], F_D_sq[to_plot, k_index], yerr=F_D_sq_unc[to_plot, k_index],
                marker='.', linestyle='none', color=color, label=plot_label)
    
    # ymins[graph_i] = np.nanmin(F_D_sq[to_plot, k_index])
    # ymaxs[graph_i] = np.nanmax(F_D_sq[to_plot, k_index])
    # assert np.isfinite(ymins[graph_i])
    # assert np.isfinite(ymaxs[graph_i])

    return_fit = False

    print('do fit', do_fit)

    if do_fit:
        # curve_fit needs the params to have similar scale (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)
        # we do not have that, so we rescale F_D by it's max value
        rescale = F_D_sq[1:, k_index].max()
        assert not np.isnan(rescale)
        F_D_sq_rescaled = F_D_sq[:, k_index] / rescale
        # F_D_sq_unc_rescaled = F_D_sq_unc[:, k_index] / rescale
        if np.isnan(F_D_sq_rescaled[1:]).sum()/F_D_sq_rescaled[1:].size == 1.0:
            print('  nan problem, skipping')
            return

        weights = np.ones_like(F_D_sq_rescaled)
        # weights = F_D_sq_unc[1:, k_index]
        # weights = F_D_sq_unc_rescaled[1:]
        weights[0] *= 1/8
        weights[1] *= 1/4
        weights[2] *= 1/2
        weights[-20:-1:2] = np.inf

        t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]), 200)

        if fit_use_flow:
            J0 = lambda x: scipy.special.j0(x)
            # func = lambda t, A, B, tau, v : A * (1 - np.exp(-t/tau) * J0(k[k_index]*v*t)) + B
            func = lambda t, A, B, D, v : A * (1 - np.exp(-t*k[k_index]**2*D) * J0(k[k_index]*v*t)) + B

            popt, pcov = scipy.optimize.curve_fit(func, t[to_plot], F_D_sq_rescaled[to_plot], sigma=weights[to_plot], p0=(F_D_sq_rescaled.max(), F_D_sq_rescaled.min(), 0.001, 0.01), maxfev=100000)#, absolute_sigma=True)
            
            D = popt[2] 
            D_unc = np.sqrt(pcov)[2][2]

            v = popt[3]
            v_unc = np.sqrt(pcov)[3][3]

            label = ''
            
            # D part of label
            if particle_diameter:
                D_SE = common.stokes_einstein_D(particle_diameter)
                label += f'$D=({common.format_val_and_unc(D/D_SE, D_unc/D_SE)})D_\mathrm{{SE}}$'
            else:
                label += f'$D={common.format_val_and_unc(D, D_unc)}\mathrm{{\mu m^2/s}}$'
            
            # v part of label
            if particle_diameter and particle_material:
                v_SE = common.stokes_einstein_v(particle_diameter, particle_material)
                if v/v_SE > 0.01:
                    label += f'\n$v=({common.format_val_and_unc(v/v_SE, v_unc/v_SE)})v_\mathrm{{SE}}$'
                else:
                    label += f'\n$v={common.format_val_and_unc(v, v_unc)}\mathrm{{\mu m/s}}$'
            else:
                label += f'\n$v={common.format_val_and_unc(v, v_unc)}\mathrm{{\mu m/s}}$'

        else: 
            func = lambda t, A, B, D: A * (1 - np.exp(-t*k[k_index]**2*D)) + B

            popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D_sq_rescaled[1:], sigma=weights[to_plot], p0=(F_D_sq_rescaled.max(), F_D_sq_rescaled.min(), 0.1), maxfev=100000)#, absolute_sigma=True)
            
            D = popt[2]
            D_unc = np.sqrt(pcov)[2][2]

            # label = f'fit $D={common.format_val_and_unc(D, D_unc)}$'
            if particle_diameter:
                D  = D / common.stokes_einstein_D(particle_diameter)
                D_unc  = D_unc / common.stokes_einstein_D(particle_diameter)
                label = f'fit $D=({D:.2g}±{D_unc:.2g})D_\mathrm{{SE}}$'
            else:
                label = f'fit $D={D:.2g}±{D_unc:.2g}\mathrm{{\mu m}}$'
        

        A = popt[0] * rescale
        B = popt[1] * rescale
        dA = np.sqrt(pcov)[0, 0] * rescale
        dB = np.sqrt(pcov)[1, 1] * rescale

        if A > 0:

            ax.plot(t_theory, func(t_theory, *popt)*rescale, color=common.FIT_COLOR, label=label, zorder=-1)

            if D_unc/D < D_ERROR_THRESH:
                return_fit = True
            else:
                print(f'  not saving D, D_unc/D={D_unc/D:.2f}')

        else:
            print('  not saving, A(q) negative')


    ax.set_title(fr'$k={k[k_index]:.2f}\mathrm{{\mu m^{{-1}}}}$ ($\approx{2*np.pi/k[k_index]:.2f}\mathrm{{\mu m}}$)')
    ax.set_xlabel('$t$ (s)')
    ax.set_ylabel('$|\Delta I(q, t)|^2$')

    ymin = F_D_sq[1:, k_index].min()/1.02
    # ymax = F_D_sq[int(F_D_sq.shape[0]*0.7), k_index] * 1.3
    ymax = F_D_sq[:, k_index].max()*1.02
    
    if particle_diameter:
        for D_mult in THEORY_PLOT_DS:
            yspan = ymax - ymin
            B = ymin + yspan/5
            A = yspan / 3
            
            D = common.stokes_einstein_D(particle_diameter) * D_mult
            
            theory = A * (1 - np.exp(-k[k_index]**2 * t_theory * D)) + B
            ax.plot(t_theory, theory, linestyle='dotted', label='$D=' f'{D_mult}' + '*D_\mathrm{Stokes-Einstein}$')

    # ax.set_ylim(ymin, ymax)


    # legend_handles, legend_labels = ax.get_legend_handles_labels()
    # if len(legend_handles):
    
        
    if SEMILOGX:
        ax.semilogx()
    # if log_y:
    #     ax.semilogy()

    if return_fit:
        return D, D_unc, k[k_index]
        

def show(file, axs, k, F_D_sq, F_D_sq_unc, t, sigma, pixel, num_displayed_ks, 
         NAME=None, channel=None, do_fit=True, k_index_offset=0, k_index_end_early=0, particle_diameter=None,
         color=None, plot_label=None, fit_use_flow=False, log_y=False, save_data=True,
         legend_fontsize=9, particle_material=None):

    real_ks = []

    Ds_for_saving = []
    D_uncs_for_saving = []
    ks_for_saving = []

    every_nth_k = (k.size - k_index_offset) // (num_displayed_ks + k_index_end_early)

    ymins = np.full(num_displayed_ks, np.nan)
    ymaxs = np.full(num_displayed_ks, np.nan)

    for graph_i in range(num_displayed_ks):
        k_index = k_index_offset + graph_i * every_nth_k

        real_ks.append(k[k_index])
        if not do_fit:
            print(f'k={k[k_index]:.2f}')

        ax = axs[graph_i]

        if (nanfrac := common.nanfrac(F_D_sq[1:, k_index])) > 0.5:
            print(f'  nanfrac = {nanfrac}')
            print('  SKIPPING')
            continue

        ret = plot_single_graph(ax, t, k, F_D_sq, F_D_sq_unc, k_index, plot_label=plot_label, color=color,
                          do_fit=do_fit, particle_diameter=particle_diameter, particle_material=particle_material,
                          fit_use_flow=fit_use_flow)
        
        if ret:
            D, D_unc, k_for_save = ret
            Ds_for_saving.append(D)
            D_uncs_for_saving.append(D_unc)
            ks_for_saving.append(k_for_save)
        
        ax.legend(fontsize=legend_fontsize)
    
    if save_data:
        common.save_data(f'visualisation/data/Ds_from_DDM_{file}',
                Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, ks=ks_for_saving,
                NAME=NAME, channel=channel)
    # common.save_data(f'DDM/data/A_B_of_q_{file}.npz',
    #          A=A_of_q, B=B_of_q, q=q,
    #          pack_frac_given=data.get('pack_frac_given'), particle_diameter=data.get('particle_diameter'))
    
    # common.save_data(f'isf/data/DDM_{file}.npz',
    #                  t=t, F=DDM_f, F_unc=DDM_f_unc, k=real_ks,
    #                  particle_diameter=data.get('particle_diameter'))

    return ymins, ymaxs, t.min(), t.max()

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'ddm_'):
        # we load this stuff here, not in `show()` as we normally do, so that we can use `show()` in the live situation
        data = common.load(f'DDM/data/ddm_{file}.npz')
        k          = data['k']
        F_D_sq     = data['F_D_sq']
        F_D_sq_unc = data['F_D_sq_unc']
        t          = data['t']

        sigma = data.get('particle_diameter')
        pixel = data['pixel_size']
        NAME  = data.get('NAME')
        channel = data.get('channel')
        
        num_displayed_ks = 14
        k_index_offset   =  0
        
        if 'psiche' in file:
            k_index_offset = k.size // 2
            num_displayed_ks = 7

        fig, axs = plt.subplots(1, num_displayed_ks, figsize=(num_displayed_ks*3.5, 4))

        show(file, axs, k, F_D_sq, F_D_sq_unc, t, sigma, pixel, NAME=NAME, channel=channel,
             num_displayed_ks=num_displayed_ks, k_index_offset=k_index_offset,
            particle_diameter=data.get('particle_diameter'), particle_material=data.get('particle_material'),
            fit_use_flow=FIT_USE_FLOW, do_fit=DO_FIT)
        
        
        suptitle = common.name(file) + ' ' + str(NAME)
        fig.suptitle(f'{suptitle}, $\sigma={sigma}$, pixel$={pixel}$')
        
        common.save_fig(fig, f'DDM/figures_png/ddm_{file}.png')


        # fig_2D, ax_2D = plt.subplots(1, 1)
        # # print(data['F_D_sq_noradial'])
        # ax_2D.imshow(data['F_D_sq_noradial'][10, :, :])
        # common.save_fig(fig_2D, f'DDM/figures_png/ddm_2d_{file}.png')