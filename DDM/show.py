import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

for file in common.files_from_argv('DDM/data', 'ddm_'):
    data = common.load(f'DDM/data/ddm_{file}.npz')
    k         = data['k']
    F_D       = data['F_D']
    t         = data['t']

    target_ks = (0.14, 1.3, 4)

    fig, axs = plt.subplots(1, len(target_ks), figsize=(len(target_ks)*3, 3.5))

    for graph_i in range(len(target_ks)):
        target_k = target_ks[graph_i]

        k_index = np.argmax(k > target_k)

        ax = axs[graph_i]

        ax.scatter(t[1:], F_D[1:, k_index], s=3)

        
        # fitting_points = common.exponential_integers(1, F_D.shape[0]-1)
        func = lambda t, A, B, tau : A * (1 - np.exp(-t/tau)) + B
        popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D[1:, k_index], p0=(1e10, 1e9, 0.1), maxfev=10000)
        t_theory = np.logspace(np.log10(t[1]), np.log10(t[-1]))
        D = 1 / (popt[2] * k[k_index]**2)
        ax.plot(t_theory, func(t_theory, *popt), color='black', label=f'fit $D={D:.3f}$\n$A=${popt[0]:.0g}, $B=${popt[1]:.0g}')


        ax.semilogx()
        # ax.semilogy()
        ax.set_title(f'$k={k[k_index]:.2f}$')
        ax.legend(fontsize=8)
    
    common.save_fig(fig, f'DDM/figures_png/ddm_{file}.png')