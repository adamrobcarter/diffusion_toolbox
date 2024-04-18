import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings
import scipy.optimize

subplot_i = 0

for type in ['DDMed']:
    for file in sys.argv[1:]:

        f_Ds_for_saving = []
        f_D_uncs_for_saving = []
        f_ks_for_saving = []
        
        Fs_Ds_for_saving = []
        Fs_D_uncs_for_saving = []
        Fs_ks_for_saving = []
        
        f_Ds_for_saving_short = []
        f_D_uncs_for_saving_short = []
        f_ks_for_saving_short = []
        
        Fs_Ds_for_saving_short = []
        Fs_D_uncs_for_saving_short = []
        Fs_ks_for_saving_short = []

        f_Ds_for_saving_long = []
        f_D_uncs_for_saving_long = []
        f_ks_for_saving_long = []
        
        Fs_Ds_for_saving_long = []
        Fs_D_uncs_for_saving_long = []
        Fs_ks_for_saving_long = []

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
        
        Fs_skold_all = F_all / F0_all
        k_skold_all = k_all / np.sqrt(F0_all)

        num_ks = k_all.shape[1]

        target_ks = (0.1, 0.14, 0.5, 1.3, 2, 4, 8)
        target_ks = list(np.logspace(np.log10(0.1), np.log10(8), 20))
        # target_ks = (0.5, 1.3, 2, 4, 8)
        # target_ks = (1.3, 6)

        # fig, (top_axes, lin_axes, D_axes) = plt.subplots(3, len(target_ks), figsize=(len(target_ks)*3, 9))

        fig, (lin_axes, D_axes) = plt.subplots(2, len(target_ks), figsize=(len(target_ks)*3, 6))

        # top_axes[0].plot(k_all[0, :], F0_all, color='black')
        # top_axes[0].semilogx()
        # top_axes[0].semilogy()
        # top_axes[0].set_title('F(0, k)')

        for graph_i in range(len(target_ks)):
            target_k = target_ks[graph_i]

            k_index = np.argmax(k_all[0, :] > target_k)
            k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

            k        = k_all      [0, k_index]
            f     = FoS_all    [:, k_index]
            f_unc = FoS_unc_all[:, k_index]
            f_unc = FoS_unc_all[:, k_index]
            F_D      = F_D_all    [:, k_index]

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
            # ax2 = top_axes[graph_i]

            label = fr"$k={k:.2f}\mathrm{{\mu m}}$ ($L\approx{2*np.pi/k:.1f}\mathrm{{\mu m}}$)"

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
            ax.errorbar(t2, Fs, yerr=Fs_unc, color='tab:orange', linestyle='', alpha=0.3)
            ax.errorbar(t,  f,  yerr=f_unc,  color='tab:blue'  , linestyle='', alpha=0.2)
            # ax2.errorbar(t , FoF0, yerr=FoF0_unc, color='tab:blue'  , linestyle='', alpha=0.2)
            # F_bad   = (3*FoF0_unc)   > FoF0
            # print('a', F_bad.shape)
            # print('b', F_bad.shape)


            
            # f_bad   = f_unc*4 > f
            # F_s_bad = Fs_unc * 4 > Fs
            f_noise   = f  < 3e-2
            F_s_noise = Fs < 1.7e-2
            f_toolong   = t > 200
            F_s_toolong = t > 400
            f_bad   = f_noise   | f_toolong
            Fs_bad = F_s_noise | F_s_toolong
            f_bad[0] = True
            Fs_bad[0] = True

            ax.scatter(t2[~Fs_bad], Fs  [~Fs_bad], label='F_s', color='tab:orange'  , s=6)
            ax.scatter(t [~f_bad  ], f[~f_bad  ], label='F/F0', color='tab:blue'    , s=6)
            # ax2.scatter(t [~F_bad  ], FoF0[~F_bad  ], label='F/F0', color='tab:blue'    , s=6)
            ax.scatter(t2[Fs_bad ], Fs  [Fs_bad ],              color='bisque'      , s=6)
            ax.scatter(t [f_bad   ], f[f_bad   ],              color='lightskyblue', s=6)
            # ax2.scatter(t [F_bad   ], FoF0[F_bad   ],              color='lightskyblue', s=6)
            # ax.scatter(t[F_s_bad ], Fs  [F_s_bad ],              color='bisque'      , s=6)
            # ax.scatter(t [~f_bad], 1-F_D[~f_bad],           color='tab:green',   s=6)
            # ax.scatter(t [ f_bad], 1-F_D[ f_bad],           color='tab:green',   s=6)
            
            
            # ax.scatter(t[:], Fs_skold[:], label=f'Fs Sköld $k={k_skold_all[0, k_index_skold]:.2f}$',            color='tab:green'      , s=6)

            # ax.set_ylim(max(Fs[~F_s_bad].max(), f[~f_bad].max())*1.1)
            # ax.set_ylim(min(Fs[~F_s_bad].min(), f[~f_bad].min())/1.1)
            
            offscreen = f <= 0
            print(f'offscreen: {offscreen.sum()/offscreen.size}')




            # fits
            # we fit to log(f(k, t)) because otherwise points near one make a much
            # larger impact to the fit than points near zero
            # print(f_unc [~f_bad].mean(), f_unc[0])
            # print(Fs_unc [~Fs_bad].mean(), Fs_unc[0])
            func = lambda t, D : np.exp(-t * k**2 * D)
            log_func = lambda t, D: np.log10( func(t, D) )
            print('av unc f ', np.log10(f_unc [~f_bad]).mean())
            print('av unc Fs', np.log10(Fs_unc [~Fs_bad]).mean())
            # log_f_unc  = np.log10(f [~f_bad ] + f_unc [~f_bad ]) - np.log10(f [~f_bad ])
            # log_Fs_unc = np.log10(Fs[~Fs_bad] + Fs_unc[~Fs_bad]) - np.log10(Fs[~Fs_bad])
            log_unc = lambda x, dx : 0.5 * np.log((x+dx)/(x-dx))
            log_f_unc  = log_unc(f , f_unc )
            log_Fs_unc = log_unc(Fs, Fs_unc)
            print('av unc f ', log_f_unc.mean())
            print('av unc Fs', log_Fs_unc.mean())
            # f_popt,  f_pcov  = scipy.optimize.curve_fit(log_func, t [~f_bad],   np.log10(f [~f_bad]),   sigma=np.log10(f_unc [~f_bad]))#,   absolute_sigma=True)
            # Fs_popt, Fs_pcov = scipy.optimize.curve_fit(log_func, t2[~Fs_bad], np.log10(Fs[~Fs_bad]), sigma=np.log10(Fs_unc[~Fs_bad]))#, absolute_sigma=True)
            f_popt,  f_pcov  = scipy.optimize.curve_fit(log_func, t [~f_bad],  np.log10(f [~f_bad]),  sigma=log_f_unc[~f_bad ],  absolute_sigma=True)
            Fs_popt, Fs_pcov = scipy.optimize.curve_fit(log_func, t2[~Fs_bad], np.log10(Fs[~Fs_bad]), sigma=log_Fs_unc[~Fs_bad], absolute_sigma=True)
            t_th = np.logspace(np.log10(t[1]), np.log10(t[-1]))
            ax.plot(t_th, func(t_th, *f_popt),  color='tab:blue', linestyle='dotted')
            ax.plot(t_th, func(t_th, *Fs_popt), color='tab:orange',   linestyle='dotted')

            f_points_short  = (~f_bad)  & (t  < 10)
            Fs_points_short = (~Fs_bad) & (t2 < 10)
            
            if f_points_short.sum() > 2:
                f_popt_short,  f_pcov_short  = scipy.optimize.curve_fit(log_func, t [f_points_short],  np.log10(f [f_points_short]),  sigma=log_f_unc [f_points_short],  absolute_sigma=True)
                ax.plot(t_th, func(t_th, *f_popt_short),  color='tab:blue', linestyle='dotted')
        
                f_Ds_for_saving_short.append(f_popt_short[0])
                f_D_uncs_for_saving_short.append(np.sqrt(f_pcov_short)[0][0])
                f_ks_for_saving_short.append(k)
            if Fs_points_short.sum() > 2:
                Fs_popt_short, Fs_pcov_short = scipy.optimize.curve_fit(log_func, t2[Fs_points_short], np.log10(Fs[Fs_points_short]), sigma=log_Fs_unc[Fs_points_short], absolute_sigma=True)
                ax.plot(t_th, func(t_th, *Fs_popt_short), color='tab:orange',   linestyle='dotted')

                Fs_Ds_for_saving_short.append(Fs_popt_short[0])
                Fs_D_uncs_for_saving_short.append(np.sqrt(Fs_pcov_short)[0][0])
                Fs_ks_for_saving_short.append(k)

            f_points_long  = (~f_bad)  & (t  > 100)
            Fs_points_long = (~Fs_bad) & (t2 > 100)
            if f_points_long.sum() > 2:
                f_popt_long,  f_pcov_long  = scipy.optimize.curve_fit(log_func, t [f_points_long,],  np.log10(f [f_points_long,]),  sigma=log_f_unc [f_points_long,],  absolute_sigma=True)
                ax.plot(t_th, func(t_th, *f_popt_long),  color='tab:blue', linestyle='dotted')
                
                f_Ds_for_saving_long.append(f_popt_long[0])
                f_D_uncs_for_saving_long.append(np.sqrt(f_pcov_long)[0][0])
                f_ks_for_saving_long.append(k)
            if Fs_points_long.sum() > 2:
                Fs_popt_long, Fs_pcov_long = scipy.optimize.curve_fit(log_func, t2[Fs_points_long,], np.log10(Fs[Fs_points_long,]), sigma=log_Fs_unc[Fs_points_long,], absolute_sigma=True)
                ax.plot(t_th, func(t_th, *Fs_popt_long), color='tab:orange',   linestyle='dotted')

                Fs_Ds_for_saving_long.append(Fs_popt_long[0])
                Fs_D_uncs_for_saving_long.append(np.sqrt(Fs_pcov_long)[0][0])
                Fs_ks_for_saving_long.append(k)


            ax.legend()
            # ax.semilogx()
            ax.semilogy()
            ax.set_title(label)
            print(~f_bad[1:5])
            print(np.argmax(f_bad[1:5]))
            print(np.argmax(f_bad[1:] ))
            ax.set_xlim(0, t[np.argmax(f_bad[1:])+5])

            ax.set_ylim(func(t[np.argmax(f_bad[1:])+5], *f_popt), 1.1)

            if np.isinf(np.sqrt(f_pcov)[0][0]):
                print(f'skipping {k:.2f}um, f_unc inf')
            else:
                print('f: D=', common.format_val_and_unc(f_popt[0], np.sqrt(f_pcov)[0][0]))
                f_Ds_for_saving.append(f_popt[0])
                f_D_uncs_for_saving.append(np.sqrt(f_pcov)[0][0])
                f_ks_for_saving.append(k)
            
            if np.isinf(np.sqrt(Fs_pcov)[0][0]):
                print(f'skipping {k:.2f}um, Fs_unc inf')
            else:
                print('Fs: D=', common.format_val_and_unc(Fs_popt[0], np.sqrt(Fs_pcov)[0][0]))
                Fs_Ds_for_saving.append(Fs_popt[0])
                Fs_D_uncs_for_saving.append(np.sqrt(Fs_pcov)[0][0])
                Fs_ks_for_saving.append(k)

            
            D_ax = D_axes[graph_i]
            
            D      = -1/(k**2 * t ) * np.log(f)
            Ds     = -1/(k**2 * t2) * np.log(Fs)
            D_unc  =  1/(k**2 * t ) / np.sqrt(f**2) * f_unc # the sqrt(**2) is needed to prevent negative errors
            Ds_unc =  1/(k**2 * t2) / np.sqrt(Fs  **2)   * Fs_unc # but remember when you do the errors properly it will be there

            D2     = -1/k**2 * np.gradient(np.log(f), t)
            D2_unc =  1/k**2 * np.gradient(1/f, t) * f_unc
            Ds2     = -1/k**2 * np.gradient(np.log(Fs), t)
            Ds2_unc =  1/k**2 * np.gradient(1/Fs, t) * Fs_unc
            
            D_ax.scatter(t [~f_bad  ], D [~f_bad  ], label='D from f', color='tab:blue'    , s=6)
            D_ax.scatter(t2[~Fs_bad], Ds[~Fs_bad], label='D from F_s' , color='tab:orange'  , s=6)
            # D_ax.scatter(t [~f_bad  ], D2[~f_bad  ], label='D2 from f' , color='tab:green'  , s=6)
            # D_ax.scatter(t2[~F_s_bad  ], Ds2[~F_s_bad  ], label='D2 from F_s' , color='tab:red'  , s=6)
            D_ax.scatter(t [ f_bad  ], D [ f_bad  ],                      color='lightskyblue', s=6)
            D_ax.scatter(t2[ Fs_bad], Ds[ Fs_bad],                      color='bisque'      , s=6)
            D_ax.errorbar(t2, Ds, yerr=Ds_unc, color='tab:orange', fmt='', alpha=0.3, linestyle='none')
            D_ax.errorbar(t , D , yerr=D_unc , color='tab:blue'  , fmt='', alpha=0.2, linestyle='none')
            
            D_ax.hlines(f_popt[0],  t.min(), t.max(), color='tab:blue', linestyle='dotted')
            D_ax.hlines(Fs_popt[0], t.min(), t.max(), color='tab:orange', linestyle='dotted')

            if f_points_short.sum() > 2:
                D_ax.hlines(f_popt_short[0],  t.min(), t.max(), color='tab:blue', linestyle='dotted')
            if Fs_points_short.sum() > 2:
                D_ax.hlines(Fs_popt_short[0], t.min(), t.max(), color='tab:orange', linestyle='dotted')

            if f_points_long.sum() > 2:
                D_ax.hlines(f_popt_long[0],  t.min(), t.max(), color='tab:blue', linestyle='dotted')
            if Fs_points_long.sum() > 2:
                D_ax.hlines(Fs_popt_long[0], t.min(), t.max(), color='tab:orange', linestyle='dotted')


            D_ax.semilogx()
            # D_ax.set_ylim(np.nanmin(D), np.nanmax(D))
            if file == 'alice0.02':
                D_ax.set_ylim(0, 0.0416*1.6)
            # if file == 'alice0.34':
            #     D_ax.set_ylim(0, 0.031*2)
            if file == 'alice0.66':
                D_ax.set_ylim(0, 0.0175*1.6)
                pass
            if file == 'eleanor0.01' or file == 'eleanor0.34':
                # D_ax.set_ylim(0, 0.08)
                # D_ax.set_ylim(0, 1)
                D_ax.set_ylim(0, np.nanmax(D))
            if file == 'eleanor0.34':
                D_ax.set_ylim(0, 0.25)
                
                # ax.relim() # tell mpl to ignore errorbars when
                # ax.autoscale_view() # calcing axis limits
                pass
            
            # D_long  = {0.34: 0.023, 0.66: 0.006}
            # D_short = {0.34: 0.033, 0.66: 0.018}
            # D_ax.axhline(y=D_short[phi], color='black', linewidth=1, linestyle=':', label=r'D from MSD')
            # D_ax.axhline(y=D_long [phi], color='black', linewidth=1, linestyle=':')
            
            D_ax.legend()
            
        plt.suptitle(fr'F or F_s (k, t), {file}')
        plt.tight_layout()

        common.save_fig(plt.gcf(), f'scattering_functions/figures_png/Fs_decay_t_{file}.png', dpi=300)
        
        common.save_data(f'visualisation/data/Ds_from_f_{file}',
                Ds=f_Ds_for_saving, D_uncs=f_D_uncs_for_saving, ks=f_ks_for_saving)
        common.save_data(f'visualisation/data/Ds_from_Fs_{file}',
                Ds=Fs_Ds_for_saving, D_uncs=Fs_D_uncs_for_saving, ks=Fs_ks_for_saving)
        common.save_data(f'visualisation/data/Ds_from_f_short_{file}',
                Ds=f_Ds_for_saving_short, D_uncs=f_D_uncs_for_saving_short, ks=f_ks_for_saving_short)
        common.save_data(f'visualisation/data/Ds_from_Fs_short_{file}',
                Ds=Fs_Ds_for_saving, D_uncs=Fs_D_uncs_for_saving, ks=Fs_ks_for_saving)
        common.save_data(f'visualisation/data/Ds_from_f_long_{file}',
                Ds=f_Ds_for_saving_long, D_uncs=f_D_uncs_for_saving_long, ks=f_ks_for_saving_long)
        common.save_data(f'visualisation/data/Ds_from_Fs_long_{file}',
                Ds=Fs_Ds_for_saving_long, D_uncs=Fs_D_uncs_for_saving_long, ks=Fs_ks_for_saving_long)