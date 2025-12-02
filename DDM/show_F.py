import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'ddm_'):
        data = common.load(f'DDM/data/ddm_{file}.npz')
        k         = data['k']
        F_D_sq    = data['F_D_sq']
        t         = data['t']

        target_ks = list(np.logspace(np.log10(0.7), np.log10(7), 5))
        # target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)
        real_ks = []

        fig, (axs, D_axs) = plt.subplots(2, len(target_ks)+1, figsize=(len(target_ks)*3.5, 6))

        D_max = 0

        for graph_i in range(len(target_ks)):
            target_k = target_ks[graph_i]

            k_index = np.argmax(k > target_k)
            real_ks.append(k[k_index])

            ax = axs[graph_i+1]

            func = lambda t, A, B, tau : A * (1 - np.exp(-t/tau)) + B
            rescale = F_D_sq[1:, k_index].max()
            # scipy doesn't really like fitting when one of the params is ~1e11 and one ~1e-1, so we rescale to sensible magnitudes for the fit
            F_D_norm = F_D_sq[:, k_index] / rescale
            if np.isnan(F_D_norm[1:]).sum()/F_D_norm[1:].size == 1.0:
                continue
            popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D_norm[1:], p0=(F_D_norm.max(), F_D_norm.min(), 0.1), maxfev=10000)
            A = popt[0] * rescale
            B = popt[1] * rescale
            t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]))
            D_POPT_INDEX = 2
            D = 1 / (popt[D_POPT_INDEX] * k[k_index]**2)
            D_unc = 1 / (k[k_index]**2 * popt[D_POPT_INDEX]**2) * np.sqrt(pcov)[D_POPT_INDEX][D_POPT_INDEX]
            # print(D, D_unc)

            F_D_sq_norm = (F_D_sq - B ) / A
            # # print('F', F_D_sq_norm[1, k_index], F_D_sq[1, k_index])
            # F_k_0 = 
            # f_k_t = 1 - F_D_sq_norm / (2 * F_k_0)
            
            
            ax.scatter(t[1:], F_D_sq_norm[1:, k_index], s=10)

            label = f'fit $D={common.format_val_and_unc(D, D_unc)}$\n$A=${popt[0]*rescale:.2g}\n$B=${popt[1]*rescale:.2g}'
            ax.plot(t_theory, np.exp(-D*k[k_index]**2*t_theory), color='black', label=label, zorder=-1)

            ax.semilogx()
            ax.semilogy()
            ax.set_ylim(1e-3, 1.1)
            ax.set_title(fr'$k={k[k_index]:.2f}$ ($\approx{2*np.pi/k[k_index]:.2f}\mathrm{{\mu m}}$)')
            ax.legend(fontsize=8)

            A = popt[0] * rescale
            B = popt[1] * rescale
            D_of_t = -1/(k[k_index]**2 * t) * np.log(1 - (F_D_sq[:, k_index] - B)/A)
            
            D_ax = D_axs[graph_i+1]
            D_ax.scatter(t[1:], D_of_t[1:], s=10)
            D_ax.semilogx()
            D_ax.hlines(D, t[1], t[-1], color='black', zorder=-1)

            print(np.nanmax(D_of_t))
            D_max = max(D_max, np.nanmax(D_of_t[np.isfinite(D_of_t)]))

        for graph_i in range(len(target_ks)):
            D_ax = D_axs[graph_i+1]
            D_ax.set_ylim(0, D_max)
        

        k_th = np.logspace(np.log10(0.7), np.log10(5), 100)
        # axs[0].plot(k_th, common.structure_factor_2d_hard_spheres(k, 0.34, 2.82))
        # axs[0].semilogx()
        # axs[0].vlines(real_ks, 0, 1.5, color='black')

        fig.suptitle(file)
        
        common.save_fig(fig, f'DDM/figures_png/ddm_{file}_F.png')