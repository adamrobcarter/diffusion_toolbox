import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings

"""

subplot_i = 0

for file in sys.argv[1:]:
    d = common.load(f"scattering_functions/data/F_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]

    d2 = common.load(f"scattering_functions/data/F_s_{file}.npz")
    t2         = d2["t"]
    Fs_all     = d2["F"]
    Fs_unc_all = d2['F_unc']
    k2_all     = d2["k"]

    assert np.array_equal(k_all, k2_all)
    print(t)
    print(t2)
    assert np.array_equal(t, t2)

    F0_all     = F_all    [0, :]
    F0_unc_all = F_unc_all[0, :]
    FoS_all = F_all / F0_all
    FoS_unc_squared = (F_unc_all / F0_all)**2 + (F_all * F0_unc_all / F0_all**2)**2
    FoS_unc_all = np.sqrt(FoS_unc_squared)

    
    Fs_skold_all = F_all / F0_all
    k_skold_all = k_all / np.sqrt(F0_all)

    num_ks = k_all.shape[1]

    target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)

    fig, (top_axes, lin_axes, D_axes) = plt.subplots(3, len(target_ks), figsize=(len(target_ks)*3, 9))

    # top_axes[0].plot(k_all[0, :], F0_all, color='black')
    # top_axes[0].semilogx()
    # top_axes[0].semilogy()
    # top_axes[0].set_title('F(0, k)')

    for graph_i in range(len(target_ks)):
        target_k = target_ks[graph_i]

        k_index = np.argmax(k_all[0, :] > target_k)
        k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

        k        = k_all      [0, k_index]
        FoF0     = FoS_all    [:, k_index]
        FoF0_unc = FoS_unc_all[:, k_index]

        Fs      = Fs_all     [:, k_index]
        Fs_unc  = Fs_unc_all [:, k_index]

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

        ax = lin_axes[graph_i]
        ax2 = top_axes[graph_i]

        label = fr"$k={k:.2f} \; (\approx{2*np.pi/k:.1f}$)"

        if np.isnan(FoF0).sum() == FoF0.size:
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
        ax.errorbar(t2, Fs  , yerr=Fs_unc,   color='tab:orange', linestyle='', alpha=0.3)
        print(t.shape, FoF0.shape)
        F_bad   = (2*FoF0_unc)   > FoF0
        # print('a', F_bad.shape)
        F_bad   = FoF0_unc   > FoF0
        # print('b', F_bad.shape)
        F_s_bad = Fs_unc > Fs

        ax.scatter(t2[~F_s_bad], Fs  [~F_s_bad], label='F_s', color='tab:orange'  , s=6)

        # ax.set_ylim(5e-4, 3)
        # ax.set_ylim(1e-4, 1e1)
        offscreen = FoF0 <= 0
        print(f'offscreen: {offscreen.sum()/offscreen.size}')

        ax.legend()
        ax.semilogx()
        ax2.semilogx()
        ax.semilogy()

        ax .set_title(label)
        ax2.set_title(label)
        ax2.set_ylim(-0.03, 0.03)

        
        D_ax = D_axes[graph_i]
        
        D      = -1/(k**2 * t ) * np.log(FoF0)
        Ds     = -1/(k**2 * t2) * np.log(Fs)
        D_unc  =  1/(k**2 * t ) / np.sqrt(FoF0**2) * FoF0_unc # the sqrt(**2) is needed to prevent negative errors
        Ds_unc =  1/(k**2 * t2) / np.sqrt(Fs  **2)   * Fs_unc # but remember when you do the errors properly it will be there
        
        D_ax.scatter(t2[~F_s_bad], Ds[~F_s_bad], label='D from F_s' , color='tab:orange'  , s=6)
        D_ax.scatter(t2[ F_s_bad], Ds[ F_s_bad],                      color='bisque'      , s=6)
        D_ax.errorbar(t2, Ds, yerr=Ds_unc, color='tab:orange', fmt='', alpha=0.3, linestyle='none')
        
        D_ax.semilogx()
        # D_ax.set_ylim(np.nanmin(D), np.nanmax(D))
        if file == 'alice0.02':
            D_ax.set_ylim(0, 0.0416*2)
        if file == 'alice0.34':
            D_ax.set_ylim(0, 0.031*2)
        if file == 'alice0.66':
            D_ax.set_ylim(0, 0.0175*2)
        
        # D_long  = {0.34: 0.023, 0.66: 0.006}
        # D_short = {0.34: 0.033, 0.66: 0.018}
        # D_ax.axhline(y=D_short[phi], color='black', linewidth=1, linestyle=':', label=r'D from MSD')
        # D_ax.axhline(y=D_long [phi], color='black', linewidth=1, linestyle=':')
        
        D_ax.legend()
        
    plt.suptitle(fr'F or F_s (k, t), {file}')
    plt.tight_layout()

    plt.savefig(f'scattering_functions/figures_png/Fs_decay_t_2_{file}.png', dpi=300)
    plt.show()