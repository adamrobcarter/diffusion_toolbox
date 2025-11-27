import countoscope_theory.nmsd
import common
import numpy as np
import warnings
import matplotlib.pyplot as plt
import countoscope_theory
import scipy.optimize

if __name__ == '__main__':
    for file in common.files_from_argv('box_counting/data/', 'counted_counts_'):
        data = common.load(f'box_counting/data/NtN0_{file}.npz')
        NtN0 = data['NtN0']
        t    = data['t']
        box_sizes = data['box_sizes']
        N    = data['N']
        N_sq = data['N_sq']
        sigma = data['particle_diameter']


        fig, ax = plt.subplots(1, 1)

        D_from_first_point = []
        L_from_first_point = []
        D_from_fit = []
        L_from_fit = []

        for box_size_index in range(NtN0.shape[0]):
            rescale_y = N_sq[box_size_index]

            ax.scatter(t[1:], NtN0[box_size_index, 1:]/rescale_y, color=common.colormap(box_size_index, 0, NtN0.shape[0]), label=fr'$L={box_sizes[box_size_index]/sigma:.2g}\sigma$')
            
            L = box_sizes[box_size_index]

            t_th = np.logspace(-0.5, 5)
            # th = N_sq[box_size_index] - N[box_size_index] * (1 - countoscope_theory.nmsd.famous_f(4*0.04*t_th/L**2)**2)
            theory = lambda t, D, N_var : N_sq[box_size_index] - N_var * (1 - countoscope_theory.nmsd.famous_f(4*D*t/L**2)**2)
            # N_var = (N_sq[box_size_index] - N[box_size_index]**2)
            th = theory(t_th, 0.04, (N_sq[box_size_index] - N[box_size_index]**2)) / rescale_y
            # ax.plot(t_th, th, color='gray', marker='none', label='sDFT (no inter)' if box_size_index==0 else None)
            # th = N_sq[box_size_index] - N[box_size_index] * (1 - countoscope_theory.nmsd.famous_f(4*0.03*t/L**2)**2)
            # th /= N_sq[box_size_index]
            # ax.plot(t[1:], th[1:], color='gray', marker='none')
            # th = N_sq[box_size_index] - N[box_size_index] * (1 - countoscope_theory.nmsd.famous_f(4*0.05*t/L**2)**2)
            # th /= N_sq[box_size_index]
            # ax.plot(t[1:], th[1:], color='gray', marker='none')

            # fit
            used_NtN0 = NtN0[box_size_index, NtN0[box_size_index, :] != 0] # idk why sometimes the last point is zero
            used_t    = t   [                NtN0[box_size_index, :] != 0]
            log_theory = lambda t, D, N_var : np.log10(theory(t, D, N_var))
            print(np.isinf(np.log10(NtN0)).sum(), 'inf')
            popt, pcov = scipy.optimize.curve_fit(log_theory, used_t, np.log10(used_NtN0), maxfev=10000, p0=(0.01, N_sq[box_size_index] - N[box_size_index]**2))
            ax.plot(t_th, theory(t_th, *popt)/rescale_y, color=common.FIT_COLOR, marker='none', label='sDFT (no inter) fit' if box_size_index==0 else None)
            D_from_fit.append(popt[0])
            L_from_fit.append(L)

            ##
            tau = 4*0.0436*t/box_sizes[box_size_index]**2

            # short = tau < 0.5
            # approx = N_sq[box_size_index] - N[box_size_index] * (2/np.sqrt(np.pi) * np.sqrt(tau[short]) - 1/np.pi * tau[short])
            # approx /= N_sq[box_size_index]
            # ax.plot(t[short], approx, color='white', marker='none')

            # short = tau < 0.5
            # approx = N_sq[box_size_index] - N[box_size_index] * (2/np.sqrt(np.pi) * np.sqrt(tau[short]))
            # approx /= N_sq[box_size_index]
            # ax.plot(t[short], approx, color='blue', marker='none')

            # long = tau > 2
            # approx = N_sq[box_size_index] - N[box_size_index] * (1 - 1/(np.pi*tau[long]) + 1/(3*np.pi*tau[long]**2))
            # approx /= N_sq[box_size_index]
            # ax.plot(t[long], approx, color='white', marker='none')

            # long = tau > 2
            # approx = N_sq[box_size_index] - N[box_size_index] * (1 - 1/(np.pi*tau[long]))
            # approx /= N_sq[box_size_index]
            # ax.plot(t[long], approx, color='blue', marker='none')

            tau_first_point = np.pi/4 * (N_sq[box_size_index] - NtN0[box_size_index, 1])**2
            # a = N[box_size_index] / np.pi
            # b = - N[box_size_index] * 2 / np.sqrt(np.pi)
            # c = N_sq[box_size_index] - NtN0[box_size_index, 1]
            # x = -b - np.sqrt(b**2 - 4*a*c) / (2*a)
            # tau_first_point = x**2
            D_first_point = tau_first_point * L**2 / (4 * t[1] * N[box_size_index]**2)
            D_from_first_point.append(D_first_point)
            L_from_first_point.append(L)

            print(fr'L={L/sigma:.1f}\sigma')
            print(f'  <N>={N[box_size_index]:.2g}, <N^2>={N_sq[box_size_index]:.2g}, <N>+<N^2>={N[box_size_index]+N[box_size_index]**2:.2g}')
            print(f'  first D={D_first_point:.4f}')

            ax.hlines(N[box_size_index]**2/rescale_y, t_th.min(), t_th.max(), color='grey', linestyle=':', zorder=-1, label=r'$\langle N \rangle^2$' if box_size_index==0 else None)
            # ax.hlines((N_sq[box_size_index]-N[box_size_index])/rescale_y, t_th.min(), t_th.max(), color='grey', linestyle='dashed', zorder=-1, label=r'$\langle N^2 \rangle - \langle N \rangle$' if box_size_index==0 else None)

        ax.set_ylabel(r'$\langle N(t)N(0) \rangle / \langle N^2 \rangle$')
        ax.set_xlabel('$t$')
        ax.semilogx()
        ax.semilogy()
        ax.legend(fontsize=8)
        
        common.save_data(f'visualisation/data/Ds_from_NtN0_first_{file}',
                Ds=D_from_first_point, D_uncs=np.zeros_like(D_from_first_point), Ls=L_from_first_point,
                particle_diameter=data.get('particle_diameter'),
                pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))
        common.save_data(f'visualisation/data/Ds_from_NtN0_fit_{file}',
                Ds=D_from_fit, D_uncs=np.zeros_like(D_from_fit), Ls=L_from_fit,
                particle_diameter=data.get('particle_diameter'),
                pixel_size=data.get('pixel_size'), window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))

        common.save_fig(fig, f'box_counting/figures_png/NtN0_{file}.png')


