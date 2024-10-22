import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings
import matplotlib.cm
import scipy.optimize

subplot_i = 0

PRESENT_SMALL = False
LOGARITHMIC_Y = False

figsize = (5, 4)
if PRESENT_SMALL:
    figsize = (3.5, 3.2)

def go(file, SHOW_FIT=False, export_destination=None):
    d = common.load(f"scattering_functions/data/F_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]

    # d2 = common.load(f"scattering_functions/data/F_s_{file}.npz")
    # t2         = d2["t"]
    # Fs_all     = d2["F"]
    # Fs_unc_all = d2['F_unc']
    # k2_all     = d2["k"]

    # assert np.array_equal(k_all, k2_all)
    # assert np.array_equal(t, t2)

    for type in ['f']:

        F0_all     = F_all    [0, :]
        F0_unc_all = F_unc_all[0, :]
        f_all = F_all / F0_all
        FoS_unc_squared = (F_unc_all / F0_all)**2 + (F_all * F0_unc_all / F0_all**2)**2
        f_unc_all = np.sqrt(FoS_unc_squared)

        
        Fs_skold_all = F_all / F0_all
        k_skold_all = k_all / np.sqrt(F0_all)

        num_ks = k_all.shape[1]

        target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)
        if PRESENT_SMALL:
            target_ks = (0.2, 0.8, 2.4)

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # top_axes[0].plot(k_all[0, :], F0_all, color='black')
        # top_axes[0].semilogx()
        # top_axes[0].semilogy()
        # top_axes[0].set_title('F(0, k)')

        for graph_i in range(len(target_ks)):
            target_k = target_ks[graph_i]

            k_index = np.argmax(k_all[0, :] > target_k)
            k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

            k     = k_all    [0, k_index]
            f     = f_all    [:, k_index]
            f_unc = f_unc_all[:, k_index]

            # Fs     = Fs_all     [:, k_index]
            # Fs_unc = Fs_unc_all [:, k_index]

            # k_skold  = k_skold_all [0, k_index_skold]
            Fs_skold = Fs_skold_all[:, k_index_skold]
            
            # top_axes[0].vlines(k, np.nanmin(F0_all), np.nanmax(F0_all), color='grey')

            # F = F_all[:, k_index]
            # F0 = F0_all[k_index]

            # print("ks", k, k_s)
            # assert np.abs(k - k_s)/k < 0.07

            # assert k_index != 0, "k_index == 0, which probably means everything is 1 because of the normalisation"
            if k_index == 0:
                warnings.warn("k_index == 0, which probably means everything is 1 because of the normalisation")
            
            # S[S < 0.02] = np.nan
            # S_s[S_s < 0.02] = np.nan

            label = fr"$k={k:.2f}\mathrm{{\mu m}}^{{-1}}$"
            # LABEL += "($L\approx{2*np.pi/k:.1f}\mathrm{{\mu m}}$)"

            if np.isnan(f).sum() == f.size:
                print(f'all nan at k={k:.1f}')
                continue

            # D, t_fit, F_fit, fit_end = individual_exponential_fit_over_t(dts, S, k)
            
            # small_decay_plot(ax, False, dts, S, t_fit, F_fit, fit_end,
            #                     title=label, show_y_labels=True)

            # data = np.load(f'msd_{phi}.npz')
            # t_msd = data['t']   * units.second
            # msd   = data['msd'] * units.micrometer**2

            # data4 = np.load(f'msd_{phi}.npz')
            # msd4  = data4['msd'] * units.micrometer**4

            # ax.plot(t_msd, np.exp(- msd * k**2 / 4), color='tab:red', label='M & T')

            # ax.plot(t_msd, 1 - 0.5 * k**2 * msd, color='tab:purple', label='taylor')

            # a2 = (msd4 - 3 * msd**2 ) / (5 * msd**2)
            # F_MnU = np.exp(- msd * k**2 / 4) * (1 + a2 * (k**2 / 4 * msd)**2/2 )
            # ax.plot(t_msd, F_MnU, color='tab:brown', label='M & U')
            # plt.legend()

            # ax.scatter(t_Fs, F_s_this, label='F_s', color='tab:orange')
            # ax.scatter(t_F , F_this  , label='F/S', color='tab:blue'  )



            f_noise   = f  < 3e-2
            # F_s_noise = Fs < 1.7e-2
            f_toolong   = t > 200
            F_s_toolong = t > 400
            f_bad   = f_noise   | f_toolong
            # F_s_bad = F_s_noise | F_s_toolong
            f_bad[0] = True
            # F_s_bad[0] = True
            

            if type == 'Fs':
                pass
            #     t_for_plot = t2
            #     F_for_plot = Fs
            #     bad_for_plot = F_s_bad
            #     F_unc_for_plot = Fs_unc
            elif type == 'f':
                t_for_plot = t
                F_for_plot = f
                bad_for_plot = f_bad
                F_unc_for_plot = f_unc

            color = common.colormap_cool(graph_i, 0, len(target_ks))
            ax.errorbar(t_for_plot, F_for_plot, yerr=F_unc_for_plot, linestyle='', alpha=0.3, color=color)
            print(t.shape, f.shape)

            ax.scatter(t_for_plot[~bad_for_plot], F_for_plot[~bad_for_plot], label='observations' if graph_i==len(target_ks)-1 else None, s=15, color=color)
            ax.scatter(t_for_plot[ bad_for_plot], F_for_plot[ bad_for_plot], alpha=0.2, s=15, color=color)
        
        
            func = lambda t, D : np.exp(-t * k**2 * D)
            log_func = lambda t, D: np.log10( func(t, D) )
            popt, pcov = scipy.optimize.curve_fit(log_func, t_for_plot[~bad_for_plot], np.log10(F_for_plot[~bad_for_plot]), sigma=np.log10(F_unc_for_plot[~bad_for_plot]))#,   absolute_sigma=True)
            t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
            theory_curve = func(t_th, popt)
            if SHOW_FIT:
                ax.plot(t_th, theory_curve, color=common.FIT_COLOR, zorder=-1, label='fit' if graph_i==len(target_ks)-1 else None)

            # ax.set_ylim(5e-4, 3)
            offscreen = f <= 0
            print(f'offscreen: {offscreen.sum()/offscreen.size}')

            if LOGARITHMIC_Y:
                t_index_for_text = int(t_th.size * (3-graph_i) / 4)
            else:
                t_index_for_text = int(350 / k_index)
            angle = np.tan(np.gradient(theory_curve, t_th)[t_index_for_text]) * 180/np.pi
            # L_label = rf'$k={k[k_index]:.1f}\mathrm{{\mu m}}^{{-1}}$'
            ax.text(t_th[t_index_for_text+0]/1.5, theory_curve[t_index_for_text+0]/1.5, label,
                    horizontalalignment='center', color=color, fontsize=9,
                    transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

        ax.legend(fontsize=7)
        ax.semilogx()
        if LOGARITHMIC_Y:
            ax.semilogy()
            ax.set_ylim(1e-2, 1.3e0)
        else:
            ax.set_ylim(-0.05, 1.05)

        if not PRESENT_SMALL:
            ax.set_title(f'$F_s(\Delta t)$, {file}')
        ax.set_xlabel('$\Delta t$ (s)')
        ylabel = '$F_s(k, \Delta t)$' if type == 'Fs' else '$f(k, \Delta t)$' if type == 'f' else None
        ax.set_ylabel(ylabel)

        fileprefix = 'fkt' if type == 'f' else 'Fs'
        filename = f'{fileprefix}_decay_overlayed_{file}'
        if SHOW_FIT:
            filename += '_fit'

        # common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/{fileprefix}_{file}.png', dpi=300, hide_metadata=True)
        if export_destination:
            common.save_fig(fig, export_destination, hide_metadata=True)
        common.save_fig(fig, f'scattering_functions/figures_png/{filename}.png', dpi=300)

        
if __name__ == '__main__':
    for file in common.files_from_argv("scattering_functions/data/", "F_"):
        go(file)