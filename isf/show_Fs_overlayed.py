import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings
import matplotlib.cm
import scipy.optimize
import visualisation.Ds_overlapped

subplot_i = 0

PRESENT_SMALL = False
LABELS_ON_PLOT = False
SHOW_THEORY = True
S_OF_K_IS_ONE = True

# PLOT_FUNCTION = lambda f: 1 - f
# PLOT_FUNCTION_ERR = lambda df: df
# PLOT_FUNCTION_NAME = '$1 - f(k, t)$'
# LOGARITHMIC_Y = True

PLOT_FUNCTION = lambda f: f
PLOT_FUNCTION_ERR = lambda df: df
PLOT_FUNCTION_NAME = '$f(k, t)$'
PLOT_FUNCTION_NAME = '$F(k, t)/F(k, 0)$'
LOGARITHMIC_Y = False


def go(file, file_i, ax, target_ks, show_fit=False, colormap=common.colormap_cool,
       show_theory=False):
    d = common.load(f"isf/data/F_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]

    # d2 = common.load(f"isf/data/F_s_{file}.npz")
    # t2         = d2["t"]
    # Fs_all     = d2["F"]
    # Fs_unc_all = d2['F_unc']
    # k2_all     = d2["k"]

    # assert np.array_equal(k_all, k2_all)
    # assert np.array_equal(t, t2)

    D_MSD, sigma, phi = visualisation.Ds_overlapped.get_D0(file)

    for f_type in ['f']:

        F0_all     = F_all    [0, :]
        F0_unc_all = F_unc_all[0, :]
        f_all = F_all / F0_all
        FoS_unc_squared = (F_unc_all / F0_all)**2 + (F_all * F0_unc_all / F0_all**2)**2
        f_unc_all = np.sqrt(FoS_unc_squared)

        
        # Fs_skold_all = F_all / F0_all
        # k_skold_all = k_all / np.sqrt(F0_all)

        # num_ks = k_all.shape[1]



        # top_axes[0].plot(k_all[0, :], F0_all, color='black')
        # top_axes[0].semilogx()
        # top_axes[0].semilogy()
        # top_axes[0].set_title('F(0, k)')

        for graph_i in range(len(target_ks)):
            target_k = list(reversed(target_ks))[graph_i] # reversed is for zorder of plots

            k_index = np.argmax(k_all > target_k)
            # k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

            k     = k_all    [k_index]
            f     = f_all    [:, k_index]
            f_unc = f_unc_all[:, k_index]

            # assert k_index != 0, "k_index == 0, which probably means everything is 1 because of the normalisation"
            if k_index == 0:
                warnings.warn("k_index == 0, which probably means everything is 1 because of the normalisation")
            
            # S[S < 0.02] = np.nan
            # S_s[S_s < 0.02] = np.nan
            # label = None
            # if file_i == 0:
            #     label = fr"$k={k:.2g}\mathrm{{\mu m}}^{{-1}}$"
            # LABEL += "($L\approx{2*np.pi/k:.1f}\mathrm{{\mu m}}$)"

            if np.isnan(f).sum() == f.size:
                print(f'all nan at k={k:.1f}')
                continue

            f_noise   = f  < 3e-2
            # F_s_noise = Fs < 1.7e-2
            f_toolong   = t > 200
            F_s_toolong = t > 400
            f_bad   = f_noise   | f_toolong
            # F_s_bad = F_s_noise | F_s_toolong
            f_bad[0] = True
            # F_s_bad[0] = True
            

            if f_type == 'Fs':
                pass
            #     t_for_plot = t2
            #     F_for_plot = Fs
            #     bad_for_plot = F_s_bad
            #     F_unc_for_plot = Fs_unc
            elif f_type == 'f':
                t_for_plot = t
                F_for_plot = f
                bad_for_plot = f_bad
                F_unc_for_plot = f_unc

            if 'longer' in file:
                print('cropping start of longer')
                thresh = 2.5e3
                F_for_plot     = F_for_plot    [t_for_plot > thresh]
                F_unc_for_plot = F_unc_for_plot[t_for_plot > thresh]
                bad_for_plot   = bad_for_plot  [t_for_plot > thresh]
                t_for_plot     = t_for_plot    [t_for_plot > thresh]

            color = colormap(graph_i/len(target_ks))
            
            if LABELS_ON_PLOT or file_i > 0:
                label = None
            else:
                
                # label = fr'$k={k:.2g}\mathrm{{\mu m^{{-1}}}}$'
                label = fr'$k\sigma={k*sigma:.2g}$'
            ax.errorbar(t_for_plot, PLOT_FUNCTION(F_for_plot), yerr=PLOT_FUNCTION_ERR(F_unc_for_plot), linestyle='', alpha=0.5, color=color)

            ax.scatter(t_for_plot[~bad_for_plot], PLOT_FUNCTION(F_for_plot[~bad_for_plot]), label=label, s=15, color=color)
            ax.scatter(t_for_plot[ bad_for_plot], PLOT_FUNCTION(F_for_plot[ bad_for_plot]), alpha=1, s=15, color=color)
        
        
            if show_fit:
                func = lambda t, D : np.exp(-t * k**2 * D)
                log_func = lambda t, D: np.log10( func(t, D) )
                popt, pcov = scipy.optimize.curve_fit(log_func, t_for_plot[~bad_for_plot], np.log10(F_for_plot[~bad_for_plot]), sigma=np.log10(F_unc_for_plot[~bad_for_plot]))#,   absolute_sigma=True)
                t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
                theory_curve = func(t_th, popt)
                label = '$\mathrm{{exp}}(-k^2D(k)t)$' if graph_i==len(target_ks)-1 else None
                ax.plot(t_th, PLOT_FUNCTION(theory_curve), color=common.FIT_COLOR, zorder=-1, label=label)

        
            if show_theory:
                if S_OF_K_IS_ONE:
                    S_of_k = 1
                    label = '$\mathrm{exp}(-k^2 t D_\mathrm{self})$' if graph_i==len(target_ks)-1 else None
                else:
                    S_of_k = common.structure_factor_2d_hard_spheres(k, phi, sigma)
                    label = '$\mathrm{exp}(-k^2 t D_\mathrm{self}/S(k))$' if graph_i==len(target_ks)-1 else None

                func = lambda t, D : np.exp(-t * k**2 * D / S_of_k)
                t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
                theory_curve = func(t_th, D_MSD)
                # theory_color = color
                theory_color = common.FIT_COLOR
                ax.plot(t_th, PLOT_FUNCTION(theory_curve), color=theory_color, zorder=-1, label=label)

            # ax.set_ylim(5e-4, 3)
            offscreen = f <= 0
            # print(f'offscreen: {offscreen.sum()/offscreen.size}')

            
            if LABELS_ON_PLOT:
                if LOGARITHMIC_Y:
                    t_index_for_text = int(t_th.size * (3-graph_i) / 4)
                else:
                    t_index_for_text = int(350 / (k_index+10))
                    
                    # this plots a straight line where we'll put the labels
                    tt = t_th
                    c = 0.35
                    m = 0.25
                    label_line = m * np.log10(tt) + c
                    # ax.plot(tt, label_line)

                    # t_index_for_text = np.argmax(F_for_plot < label_line)
                    t_index_for_text = np.argmax(theory_curve < label_line)

                # t_index_for_text = min(t_index_for_text, t_th.size-1)

                # t_index_for_text = np.argmax(f<0.5)

                angle = np.arctan(np.gradient(theory_curve, t_th)[t_index_for_text]) * 180/np.pi
                # L_label = rf'$k={k[k_index]:.1f}\mathrm{{\mu m}}^{{-1}}$'
                ax.text(t_th[t_index_for_text]*1.1, theory_curve[t_index_for_text]+0.03, label,
                        horizontalalignment='center', color=color, fontsize=9,
                        transform_rotates_text=True, rotation=angle, rotation_mode='anchor')

            else:
                ax.legend(fontsize=7)

        ax.semilogx()
        if LOGARITHMIC_Y:
            ax.semilogy()
            ax.set_ylim(3e-4, 1.3e0)
        else:
            ax.set_ylim(-0.05, 1.05)

        ax.set_xlabel('$t$ (s)')
        print('type', f_type)
        # ylabel = '$F_s(k, t)$' if f_type == 'Fs' else '$f(k, t)$' if f_type == 'f' else None
        ylabel = PLOT_FUNCTION_NAME
        ax.set_ylabel(ylabel)



figsize = (5, 4)
if PRESENT_SMALL:
    figsize = (3.5, 3.2)
fig, ax = plt.subplots(1, 1, figsize=figsize)
        
        
if __name__ == '__main__':
    for file_i, file in enumerate(common.files_from_argv("isf/data/", "F_")):
        
        target_ks = np.logspace(np.log10(0.002), np.log10(2), 10)
        target_ks = np.logspace(np.log10(0.02), np.log10(4), 10)
        # target_ks[-3] *= 1.5
        #0.001, 0.14, 0.5, 1.3, 2, 4, 8)
        if PRESENT_SMALL:
            target_ks = (0.2, 0.8, 2.4, 4.8)

        # if 'longer' in file:
        #     target_ks = (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

        go(file, file_i, ax, target_ks=target_ks,
           show_theory=SHOW_THEORY)

        if not PRESENT_SMALL:
            ax.set_title(f'{PLOT_FUNCTION_NAME}, {file}')

    filename = f'f_decay_overlayed_{file}'

    if S_OF_K_IS_ONE:
        filename += '_Sk_is_one'
    common.save_fig(fig, f'isf/figures_png/{filename}.png', dpi=300)