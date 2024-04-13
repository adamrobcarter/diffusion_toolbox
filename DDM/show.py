import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

for file in common.files_from_argv('DDM/data', 'ddm_'):
    data = common.load(f'DDM/data/ddm_{file}.npz')
    k         = data['k']
    F_D_sq    = data['F_D_sq']
    t         = data['t']

    target_ks = list(np.logspace(np.log10(0.25), np.log10(7), 8))
    # target_ks = list(np.logspace(np.log10(0.4), np.log10(7), 7))
    # target_ks = (0.28, 0.38, 0.5, 1.3, 2, 4, 8)
    # target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)
    real_ks = []

    fig, (axs, D_axs) = plt.subplots(2, len(target_ks)+1, figsize=(len(target_ks)*3.5, 6))

    D_max = 0

    Ds_for_saving = []
    D_uncs_for_saving = []
    ks_for_saving = []

    for graph_i in range(len(target_ks)):
        target_k = target_ks[graph_i]

        k_index = np.argmax(k > target_k)
        real_ks.append(k[k_index])

        ax = axs[graph_i+1]

        ax.scatter(t[1:], F_D_sq[1:, k_index], s=10)


        func = lambda t, A, B, tau : A * (1 - np.exp(-t/tau)) + B
        rescale = F_D_sq[1:, k_index].max()
        F_D_sq_rescaled = F_D_sq[:, k_index] / rescale
        if np.isnan(F_D_sq_rescaled[1:]).sum()/F_D_sq_rescaled[1:].size == 1.0:
            continue

        weights = np.ones_like(F_D_sq_rescaled[1:])
        weights[0] = 1/8
        weights[1] = 1/4
        weights[2] = 1/2
        # popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D_norm[1:], p0=(F_D_norm.max(), F_D_norm.min(), 0.1), maxfev=10000)
        popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D_sq_rescaled[1:], sigma=weights, p0=(F_D_sq_rescaled.max(), F_D_sq_rescaled.min(), 0.1), maxfev=10000)
        t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]))
        D_POPT_INDEX = 2
        D = 1 / (popt[D_POPT_INDEX] * k[k_index]**2)
        D_unc = 1 / (k[k_index]**2 * popt[D_POPT_INDEX]**2) * np.sqrt(pcov)[D_POPT_INDEX][D_POPT_INDEX]
        # print(D, D_unc)
        label = f'fit $D={common.format_val_and_unc(D, D_unc)}$\n$A=${popt[0]*rescale:.2g}\n$B=${popt[1]*rescale:.2g}'
        ax.plot(t_theory, func(t_theory, *popt)*rescale, color='black', label=label, zorder=-1)
        # ax.plot(t_theory, func(t_theory, popt[0], popt[1], 1/(0.03*k[k_index]**2))*rescale, color='grey', label=label, zorder=-1)

        Ds_for_saving.append(D)
        D_uncs_for_saving.append(D_unc)
        ks_for_saving.append(k[k_index])

        ax.semilogx()
        # ax.semilogy()
        ax.set_title(fr'$k={k[k_index]:.2f}$ ($\approx{2*np.pi/k[k_index]:.2f}\mathrm{{\mu m}}$)')
        ax.legend(fontsize=8)

        A = popt[0] * rescale
        B = popt[1] * rescale
        D_of_t = -1/(k[k_index]**2 * t) * np.log(1 - (F_D_sq[:, k_index] - B)/A)
        
        D_ax = D_axs[graph_i+1]
        D_ax.scatter(t[1:], D_of_t[1:], s=10)
        D_ax.semilogx()
        D_ax.hlines(D, t[1], t[-1], color='black', zorder=-1)

        D_max = max(D_max, np.nanmax(D_of_t[np.isfinite(D_of_t)]))

    for graph_i in range(len(target_ks)):
        D_ax = D_axs[graph_i+1]
        D_ax.set_ylim(0, D_max*1.1)
        # D_ax.set_ylim(0, D_max/5)

    k_th = np.logspace(np.log10(0.05), np.log10(10), 100)
    axs[0].plot(k_th, common.structure_factor_2d_hard_spheres(k_th, 0.34, 2.82))
    axs[0].semilogx()
    axs[0].vlines(real_ks, 0, 1.5, color='black')

    sigma = data['particle_diameter']
    pixel = data['pixel_size']
    fig.suptitle(f'{file}, $\sigma={sigma}$, pixel$={pixel}$')
    
    common.save_fig(fig, f'DDM/figures_png/ddm_{file}.png')
    np.savez(f'visualisation/data/Ds_from_DDM_{file}',
             Ds=Ds_for_saving, D_uncs=D_uncs_for_saving, ks=ks_for_saving)