import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import matplotlib.cm

for file in common.files_from_argv('DDM/data', 'ddm_'):
    data = common.load(f'DDM/data/ddm_{file}.npz')
    k         = data['k']
    F_D_sq    = data['F_D_sq']
    t         = data['t']

    real_ks = []

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.4))
    ax.semilogx()
    ax.set_xlabel('$\Delta t$ (s)')
    ax.set_ylabel('$d(k, \Delta t)$')

    target_ks = (0.2, 0.8, 2.4)

    for graph_i in range(len(target_ks)):
        target_k = target_ks[graph_i]

        k_index = np.argmax(k > target_k)

        # label = fr"$k={k:.2f}\mathrm{{\mu m}}^{{-1}}$"
        color = matplotlib.cm.afmhot((graph_i+0.8)/(len(target_ks)+2))
        ax.scatter(t[1:], F_D_sq[1:, k_index], s=15, color=color)

    common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/ddm_overlapped_nofit_{file}.pdf', hide_metadata=True)
    

    for graph_i in range(len(target_ks)):
        target_k = target_ks[graph_i]
        k_index = np.argmax(k > target_k)

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
        theory_curve = func(t_theory, *popt)*rescale
        ax.plot(t_theory, theory_curve, color='black', zorder=-1)
        # ax.plot(t_theory, func(t_theory, popt[0], popt[1], 1/(0.03*k[k_index]**2))*rescale, color='grey', label=label, zorder=-1)

        t_index_for_text = int(t_theory.size * (3-graph_i) / 4)
        print(np.gradient(theory_curve, t_theory)[t_index_for_text])
        angle = np.arctan(np.gradient(theory_curve, t_theory)[t_index_for_text]) * 180/np.pi
        # angle = 0
        print(angle)
        # angle = 89
        # angle += 180
        angle = 58
        # print(angle)
        L_label = rf'$k={k[k_index]:.1f}\mathrm{{\mu m}}^{{-1}}$'
        ax.text(t_theory[t_index_for_text+0]*0.7, theory_curve[t_index_for_text+0]*1, L_label,
                horizontalalignment='center', color=color, fontsize=9,
                # transform_rotates_text=True, rotation=angle, rotation_mode='anchor')
                rotation=angle, rotation_mode='anchor')

    # ax.semilogy()
    # ax.set_title(fr'$k={k[k_index]:.2f}$ ($\approx{2*np.pi/k[k_index]:.2f}\mathrm{{\mu m}}$)')
    # ax.legend(fontsize=8)

    common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/ddm_overlapped_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'DDM/figures_png/ddm_overlapped_{file}.png', dpi=200)
    
