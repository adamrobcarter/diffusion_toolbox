import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

for file in common.files_from_argv('DDM/data', 'ddm_'):
    data = common.load(f'DDM/data/ddm_{file}.npz')
    k         = data['k']
    F_D_sq    = data['F_D_sq']
    t         = data['t']

    real_ks = []

    fig, ax = plt.subplots(1, 1, figsize=(3.3, 3))


    target_k = 1

    k_index = np.argmax(k > target_k)

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


    ax.semilogx()
    # ax.semilogy()
    ax.set_title(fr'$k={k[k_index]:.2f}$ ($\approx{2*np.pi/k[k_index]:.2f}\mathrm{{\mu m}}$)')
    ax.legend(fontsize=8)

    ax.set_xlabel('$t$')
    ax.set_ylabel('$D(k, t)$')

    common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/ddm_{file}.png', hide_metadata=True)

# for graph_i in range(len(target_ks)):

#     k_th = np.logspace(np.log10(0.05), np.log10(10), 100)
#     axs[0].plot(k_th, common.structure_factor_2d_hard_spheres(k_th, 0.34, 2.82))
#     axs[0].semilogx()
#     axs[0].vlines(real_ks, 0, 1.5, color='black')

#     sigma = data['particle_diameter']
#     pixel = data['pixel_size']
#     fig.suptitle(f'{file}, $\sigma={sigma}$, pixel$={pixel}$')
    
