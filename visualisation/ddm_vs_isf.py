import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings
import scipy.optimize

subplot_i = 0

for file in sys.argv[1:]:
    d = common.load(f"scattering_functions/data/F_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]

    F0_all     = F_all    [0, :]
    F0_unc_all = F_unc_all[0, :]
    FoS_all = F_all / F0_all
    FoS_unc_squared = (F_unc_all / F0_all)**2 + (F_all * F0_unc_all / F0_all**2)**2
    FoS_unc_all = np.sqrt(FoS_unc_squared)

    F_D_all = 2 * F0_all - 2 * F_all

    d2 = common.load(f"scattering_functions/data/F_s_{file}.npz")
    t2         = d2["t"]
    Fs_all     = d2["F"]
    Fs_unc_all = d2['F_unc']
    k2_all     = d2["k"]

    assert np.array_equal(k_all, k2_all)
    assert np.array_equal(t, t2)

    d3 = common.load(f'DDM/data/ddm_{file}.npz')
    k_DDM_all  = d3['k']
    F_D_sq_all = d3['F_D_sq']
    t_DDM      = d3['t']
    
    Fs_skold_all = F_all / F0_all
    k_skold_all = k_all / np.sqrt(F0_all)

    num_ks = k_all.shape[1]

    target_ks = (0.14, 0.5, 1.3, 2, 4, 8)
    # target_ks = (0.1, 0.2, 0.28, 0.38, 0.5, 1.3, 2, 4, 8)
    target_ks = list(np.logspace(np.log10(0.05), np.log10(7), 9))
    # target_ks = (0.5, 1.3, 2, 4, 8)
    # target_ks = (1.3, 6)

    fig, (lin_axes, D_axes) = plt.subplots(2, len(target_ks), figsize=(len(target_ks)*3, 6.5))


    for graph_i in range(len(target_ks)):
        target_k = target_ks[graph_i]

        print(k_DDM_all.max())

        k_index       = np.argmax(k_all      [0, :] > target_k)
        k_index_DDM   = np.argmax(k_DDM_all  [:]    > target_k)
        k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

        k        = k_all      [0, k_index]
        FoF0     = FoS_all    [:, k_index]
        FoF0_unc = FoS_unc_all[:, k_index]
        FoF0_unc = FoS_unc_all[:, k_index]
        two_Fk0_minus_Fkt      = F_D_all    [:, k_index]

        Fs      = Fs_all     [:, k_index]
        Fs_unc  = Fs_unc_all [:, k_index]

        # k_skold  = k_skold_all [0, k_index_skold]
        Fs_skold = Fs_skold_all[:, k_index_skold]

        k_DDM  = k_DDM_all [k_index_DDM]
        F_D_sq = F_D_sq_all[:, k_index_DDM]

        # F = F_all[:, k_index]
        # F0 = F0_all[k_index]

        # assert k_index != 0, "k_index == 0, which probably means everything is 1 because of the normalisation"
        if k_index == 0:
            warnings.warn("k_index == 0, which probably means everything is 1 because of the normalisation")

        ax = lin_axes[graph_i]

        label = fr"$k={k:.2f}\mathrm{{\mu m}}$ ($L\approx{2*np.pi/k:.1f}\mathrm{{\mu m}}$)" "\n" fr" $k_\mathrm{{DDM}}={k_DDM:.2f}\mathrm{{\mu m}}$"

        if np.isnan(FoF0).sum() == FoF0.size:
            print(f'all nan at k={k:.1f}')
            continue

        # ax.scatter(t_Fs, F_s_this, label='F_s', color='tab:orange')
        # ax.scatter(t_F , F_this  , label='F/S', color='tab:blue'  )
        # ax.errorbar(t2, Fs  , yerr=Fs_unc,   color='tab:orange', linestyle='', alpha=0.3)
        # print(t.shape, FoF0.shape)
        # ax .errorbar(t , FoF0, yerr=FoF0_unc, color='tab:blue'  , linestyle='', alpha=0.2)
        # ax2.errorbar(t , FoF0, yerr=FoF0_unc, color='tab:blue'  , linestyle='', alpha=0.2)
        # F_bad   = (3*FoF0_unc)   > FoF0
        # print('a', F_bad.shape)
        # print('b', F_bad.shape)
        
        F_bad   = FoF0_unc*4 > FoF0
        F_s_bad = Fs_unc * 4 > Fs
        F_bad = FoF0 < 2e-2
        F_s_bad = Fs < 1.7e-2

        
        weights = np.ones_like(F_D_sq[1:])
        weights[0] = 1/8
        weights[1] = 1/4
        weights[2] = 1/2
        func = lambda t, A, B, tau : A * (1 - np.exp(-t/tau)) + B
        rescale = F_D_sq[1:].max()
        # scipy doesn't really like fitting when one of the params is ~1e11 and one ~1e-1, so we rescale to sensible magnitudes for the fit
        F_D_sq_norm = F_D_sq / rescale
        if np.isnan(F_D_sq_norm[1:]).sum()/F_D_sq_norm[1:].size == 1.0:
            continue
        popt, pcov = scipy.optimize.curve_fit(func, t[1:], F_D_sq_norm[1:], sigma=weights, p0=(F_D_sq_norm.max(), F_D_sq_norm.min(), 0.1), maxfev=10000)
        A = popt[0] * rescale
        B = popt[1] * rescale

        F_D_sq_norm = (F_D_sq - B ) / A

        # ax.scatter(t2[~F_s_bad], Fs  [~F_s_bad], label='F_s', color='tab:orange'  , s=6)
        # ax.scatter(t [~F_bad  ], FoF0[~F_bad  ], label='F/F0', color='tab:blue'    , s=6)
        # # ax2.scatter(t [~F_bad  ], FoF0[~F_bad  ], label='F/F0', color='tab:blue'    , s=6)
        # ax.scatter(t2[F_s_bad ], Fs  [F_s_bad ],              color='bisque'      , s=6)
        # ax.scatter(t [F_bad   ], FoF0[F_bad   ],              color='lightskyblue', s=6)
        # ax2.scatter(t [F_bad   ], FoF0[F_bad   ],              color='lightskyblue', s=6)
        # ax.scatter(t[F_s_bad ], Fs  [F_s_bad ],              color='bisque'      , s=6)
        # ax.scatter(t [~F_bad], 1-F_D[~F_bad],           color='tab:green',   s=6)
        # ax.scatter(t [ F_bad], 1-F_D[ F_bad],           color='tab:green',   s=6)
        # ax.scatter(t, two_Fk0_minus_Fkt/2, s=6, label='$F(k, 0) - F(k, t)$')
        # ax.scatter(t, two_Fk0_minus_Fkt, s=6, label='$2F(k, 0) - 2F(k, t)$')
        ax.scatter(t, Fs,   s=6, label='$F_s(k, t)$')
        ax.scatter(t, FoF0, s=6, label='$f(k, t)$')
        ax.scatter(t, 1-F_D_sq_norm, s=6, label=r'$1 - D_\mathrm{norm}(k, t)$')
        
        
        # ax.scatter(t[:], Fs_skold[:], label=f'Fs Sköld $k={k_skold_all[0, k_index_skold]:.2f}$',            color='tab:green'      , s=6)

        # ax.set_ylim(5e-4, 3)
        # ax.set_ylim(1e-4, 1e1)
        offscreen = FoF0 <= 0
        print(f'offscreen: {offscreen.sum()/offscreen.size}')

        ax.legend()
        ax.semilogx()
        # ax2.semilogx()
        ax.semilogy()

        ax .set_title(label)
        # ax2.set_title(label)
        # ax2.set_ylim(-0.03, 0.03)

        
        D_ax = D_axes[graph_i]
        
        D      = -1/(k**2 * t ) * np.log(FoF0)
        Ds     = -1/(k**2 * t2) * np.log(Fs)
        D_DDM  = -1/(k_DDM**2 * t_DDM) * np.log(1-F_D_sq_norm)
        D_unc  =  1/(k**2 * t ) / np.sqrt(FoF0**2) * FoF0_unc # the sqrt(**2) is needed to prevent negative errors
        Ds_unc =  1/(k**2 * t2) / np.sqrt(Fs  **2)   * Fs_unc # but remember when you do the errors properly it will be there

        D2     = -1/k**2 * np.gradient(np.log(FoF0), t)
        D2_unc =  1/k**2 * np.gradient(1/FoF0, t) * FoF0_unc
        
        D_ax.scatter(t2,    Ds,    label='D from $F_s$', s=6)
        D_ax.scatter(t,     D,     label='D from $f$', s=6)
        D_ax.scatter(t_DDM, D_DDM, label=r'D from $1-D_\mathrm{norm}$', s=6)
        # # D_ax.scatter(t [~F_bad  ], D2[~F_bad  ], label='D2' , color='tab:green'  , s=6)
        # D_ax.scatter(t [ F_bad  ], D [ F_bad  ],                      color='lightskyblue', s=6)
        # D_ax.scatter(t2[ F_s_bad], Ds[ F_s_bad],                      color='bisque'      , s=6)
        # D_ax.errorbar(t2, Ds, yerr=Ds_unc, color='tab:orange', fmt='', alpha=0.3, linestyle='none')
        # D_ax.errorbar(t , D , yerr=D_unc , color='tab:blue'  , fmt='', alpha=0.2, linestyle='none')
        
        D_ax.semilogx()
        # D_ax.set_ylim(np.nanmin(D), np.nanmax(D))
        if file == 'alice0.02':
            D_ax.set_ylim(0, 0.0416*1.6)
        # if file == 'alice0.34':
        #     D_ax.set_ylim(0, 0.031*2)
        if file == 'alice0.66':
            D_ax.set_ylim(0, 0.0175*1.6)
            pass
        if file == 'eleanor0.01':
            # D_ax.set_ylim(0, 0.08)
            pass
        if file == 'eleanor0.34':
            # D_ax.set_ylim(0, 0.5)
            pass
        
        # D_long  = {0.34: 0.023, 0.66: 0.006}
        # D_short = {0.34: 0.033, 0.66: 0.018}
        # D_ax.axhline(y=D_short[phi], color='black', linewidth=1, linestyle=':', label=r'D from MSD')
        # D_ax.axhline(y=D_long [phi], color='black', linewidth=1, linestyle=':')
        
        D_ax.legend()
        
    plt.suptitle(fr'F or F_s (k, t), {file}')

    common.save_fig(fig, f'visualisation/figures_png/ddm_vs_isf_{file}.png', dpi=300)