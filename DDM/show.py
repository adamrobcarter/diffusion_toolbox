import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.special
import math

FIT_USE_FLOW = False
D_ERROR_THRESH = 0.5 # D_unc/D must be smaller than this to save

def show(file, axs, k, F_D_sq, F_D_sq_unc, t, sigma, pixel, num_displayed_ks, NAME=None, channel=None, live=False):

    real_ks = []

    D_max = 0

    Ds_for_saving = []
    D_uncs_for_saving = []
    ks_for_saving = []

    A_of_q = []
    B_of_q = []
    q = []

    every_nth_k = k.size // num_displayed_ks

    for graph_i in range(num_displayed_ks):
        k_index = graph_i * every_nth_k

        real_ks.append(k[k_index])
        if not live:
            print(f'k={k[k_index]:.2f}')

        ax = axs[graph_i]

        if (nanfrac := common.nanfrac(F_D_sq[1:, k_index])) > 0.5:
            print(f'  nanfrac = {nanfrac}')
            print('  SKIPPING')
            continue

        to_plot = np.full_like(F_D_sq[:, k_index], True, dtype='bool')
        to_plot[0] = False
        anomalous = F_D_sq_unc[:, k_index] > F_D_sq[:, k_index]
        # to_plot[anomalous] = False


        if num_removed := (to_plot==False).sum()-1:
            print(f'  removed {num_removed}')

        ax.errorbar(t[to_plot], F_D_sq[to_plot, k_index], yerr=F_D_sq_unc[to_plot, k_index], marker='.', linestyle='none')

        if not live:
            # curve_fit needs the params to have similar scale (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)
            # we do not have that, so we rescale F_D by it's max value
            rescale = F_D_sq[1:, k_index].max()
            assert not np.isnan(rescale)
            F_D_sq_rescaled = F_D_sq[:, k_index] / rescale
            # F_D_sq_unc_rescaled = F_D_sq_unc[:, k_index] / rescale
            if np.isnan(F_D_sq_rescaled[1:]).sum()/F_D_sq_rescaled[1:].size == 1.0:
                print('  nan problem, skipping')
                continue

            weights = np.ones_like(F_D_sq_rescaled)
            # weights = F_D_sq_unc[1:, k_index]
            # weights = F_D_sq_unc_rescaled[1:]
            weights[0] *= 1/8
            weights[1] *= 1/4
            weights[2] *= 1/2
            weights[-20:-1:2] = np.inf

            t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]), 200)

            if FIT_USE_FLOW:
                J0 = lambda x: scipy.special.j0(x)
                # func = lambda t, A, B, tau, v : A * (1 - np.exp(-t/tau) * J0(k[k_index]*v*t)) + B
                func = lambda t, A, B, D, v : A * (1 - np.exp(-t*k[k_index]**2*D) * J0(k[k_index]*v*t)) + B

                popt, pcov = scipy.optimize.curve_fit(func, t[to_plot], F_D_sq_rescaled[to_plot], sigma=weights[to_plot], p0=(F_D_sq_rescaled.max(), F_D_sq_rescaled.min(), 0.001, 0.01), maxfev=100000)#, absolute_sigma=True)
                
                D = popt[2] 
                D_unc = np.sqrt(pcov)[2][2]

                v = popt[3]
                v_unc = np.sqrt(pcov)[3][3]
                
                label = f'$D={common.format_val_and_unc(D, D_unc)}$\n$v={common.format_val_and_unc(v, v_unc)}$'
                # label = f'fit $D={D:.5f}$, $v={v:.5f}$'

            else: 
                func = lambda t, A, B, D: A * (1 - np.exp(-t*k[k_index]**2*D)) + B

                popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D_sq_rescaled[1:], sigma=weights[to_plot], p0=(F_D_sq_rescaled.max(), F_D_sq_rescaled.min(), 0.1), maxfev=100000)#, absolute_sigma=True)
                
                D = popt[2] 
                D_unc = np.sqrt(pcov)[2][2]

                # label = f'fit $D={common.format_val_and_unc(D, D_unc)}$'
                label = f'fit $D={D:.2g}±{D_unc:.2g}$'
            

            A = popt[0] * rescale
            B = popt[1] * rescale
            dA = np.sqrt(pcov)[0, 0] * rescale
            dB = np.sqrt(pcov)[1, 1] * rescale
            A_of_q.append(A)
            B_of_q.append(B)
            q.append(k[k_index])

            if A > 0:

                ax.plot(t_theory, func(t_theory, *popt)*rescale, color=common.FIT_COLOR, label=label, zorder=-1)

                if D_unc/D < D_ERROR_THRESH:
                    Ds_for_saving.append(D)
                    D_uncs_for_saving.append(D_unc)
                    ks_for_saving.append(k[k_index])
                else:
                    print(f'  not saving D, D_unc/D={D_unc/D:.2f}')

            else:
                print('  not saving, A(q) negative')


        ax.semilogx()
        # ax.semilogy()
        ax.set_title(fr'$k={k[k_index]:.2f}$ ($\approx{2*np.pi/k[k_index]:.2f}\mathrm{{\mu m}}$)')
        ax.set_xlabel('$t$ (s)')
        ax.set_ylabel('$|\Delta I(q, t)|^2$')
        
        legend_handles, legend_labels = ax.get_legend_handles_labels()
        if len(legend_handles):
            ax.legend(fontsize=9)


    fig.suptitle(f'{common.name(file)}, $\sigma={sigma}$, pixel$={pixel}$')
    

    common.save_data(f'visualisation/data/Ds_from_DDM_{file}',
             Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, ks=ks_for_saving,
             NAME=NAME, channel=channel)
    # common.save_data(f'DDM/data/A_B_of_q_{file}.npz',
    #          A=A_of_q, B=B_of_q, q=q,
    #          pack_frac_given=data.get('pack_frac_given'), particle_diameter=data.get('particle_diameter'))
    
    # common.save_data(f'isf/data/DDM_{file}.npz',
    #                  t=t, F=DDM_f, F_unc=DDM_f_unc, k=real_ks,
    #                  particle_diameter=data.get('particle_diameter'))

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'ddm_'):
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
        fig, axs = plt.subplots(1, num_displayed_ks, figsize=(num_displayed_ks*3.5, 4))

        show(file, axs, k, F_D_sq, F_D_sq_unc, t, sigma, pixel, NAME=NAME, channel=channel,
             num_displayed_ks=num_displayed_ks)
        
        common.save_fig(fig, f'DDM/figures_png/ddm_{file}.png')