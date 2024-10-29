import numpy as np
import common
import matplotlib.pyplot as plt
import sys
import warnings, math
import scipy.optimize, scipy.signal
import matplotlib.cm

subplot_i = 0

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def show_single_F_type(file, fig, ax, num_displayed_ks, mult=False):
    # lin_short_axes, lin_axes, log_axes, D_axes, extra_axes = axes


    d = common.load(f"isf/data/F_{file}.npz")
    t         = d["t"]
    F_all     = d["F"]
    F_unc_all = d['F_unc']
    k_all     = d["k"]
    particle_diameter     = d['particle_diameter']

    every_nth_k = int(math.ceil(k_all.shape[1] / num_displayed_ks))
    every_nth_k = max(every_nth_k, 1)

    graph_i = 0

    for k_index in range(k_all.shape[1]):
        # target_k = target_ks[graph_i]
            # k_index_skold = np.argmax(k_skold_all[0, :] > target_k)

        ks = k_all[0, :]

        # k_index = np.argmax(ks >= target_k)
        k = ks[k_index]


        # print(f'k: target {target_k:.3f}, real {k:.3f}, index {k_index}, 2pi/k={2*np.pi/k:.1f}um')
        print(f'k {k:.3f}, index {k_index}, 2pi/k={2*np.pi/k:.1f}um')

        f     = F_all    [:, k_index]
        f_unc = F_unc_all[:, k_index]
        
        f /= F_all[0, k_index]
        f_unc_sq_all = (F_unc_all / F_all[0, :])**2 + (F_all * F_unc_all[0, :] / F_all[0, :]**2)**2
        # assert f_unc_sq.shape == f_unc.shape, f'{f_unc_sq.shape} {f_unc.shape}'
        f_unc = np.sqrt(f_unc_sq_all)[:, k_index]


        display = False
        if k_index % every_nth_k == 0:
            display = True

            graph_i += 1

        label = fr"$k={k:.3f}\mathrm{{\mu m}}$ ($L\approx{2*np.pi/k:.1f}\mathrm{{\mu m}}$, $k\sigma={k*particle_diameter:.2f}$)"

        if common.nanfrac(f) == 1:
            print(f'all nan at k={k:.1f} (i={graph_i})')
            continue


        noise_thresh = 1e-2 # for eleanorlong!!
        time_thresh = 200

        f_noise   = f < noise_thresh
        f_toolong = t > time_thresh
            # f_bad   = f_noise | f_toolong # f_toolong should depend on k!!!!
        f_bad   = f_noise
        f_bad[0] = True

        # new noise identification idea
        # first noise point is first point where the gradient is no longer getting more negative
        grad = np.gradient(np.log10(f), np.log10(t))
        peaks, props = scipy.signal.find_peaks(-grad, prominence=0.01)
        # prominance filters out peaks that are just noise. A smoothing filter would probably be better though
        f_bad = np.full(f.shape, False)
        
        if len(peaks):
            f_bad[peaks[0]:] = True
        else:
            print('  no peaks in -grad found?!')


        w = -np.log(f) / k**2
        w_unc = np.abs( 1/(k**2 * f) * f_unc )
            
        print('uncs', np.mean(f_unc/f), np.nanmean(w_unc/w))

        f_bad[f<0] = True

        if display:
            color = matplotlib.cm.afmhot(np.interp(graph_i, (0, num_displayed_ks), (0.2, 0.8)))
            
            ax.plot(t[~f_bad], w[~f_bad],                         color=color, linestyle='', label=f'k={k:.2f}', marker='.')
            ax.errorbar(t[~f_bad], w[~f_bad], yerr=w_unc[~f_bad], color=color, linestyle='', alpha=0.2,          marker='none')
            ax.errorbar(t[ f_bad], w[ f_bad], yerr=w_unc[ f_bad], color=color, linestyle='', alpha=0.2,          marker='.')

    ax.loglog()
    ax.legend(fontsize=6)

if __name__ == '__main__':
    for file in sys.argv[1:]:

        
        num_displayed_ks = 20
        fig, axes = plt.subplots(1, 1)
        
        show_single_F_type(file, fig, axes, num_displayed_ks)

        common.save_fig(fig, f'isf/figures_png/f_width_{file}.png', dpi=300)
        