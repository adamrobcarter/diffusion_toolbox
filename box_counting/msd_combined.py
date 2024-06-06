import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import common
import scipy.integrate
import sys
import matplotlib.cm

fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

titles = []

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors += colors # concat in case we have not enough colors, lazy
color_index = 0
if sys.argv[1] == 'alice0.02_overlapped3' and sys.argv[2] == 'alice0.02_overlapped':
    color_index += 1
if sys.argv[1] == 'alice0.02_overlapped3' and sys.argv[2] == 'alice0.02_overlapped_neg':
    color_index += 1

collapse_x = True
collapse_y = True
collapse_x = False
collapse_y = False

Ds = {}
D_uncs = {}

for file in (files := common.files_from_argv('box_counting/data', 'counted_')):
    # ax.set_prop_cycle(None) # reset colour cycle
    color = colors[color_index]
    color_index += 1

    D0_from_fits     = [{}, {}]
    D0_unc_from_fits = [{}, {}]

    LOWTIME_FIT_END = 20

    boxes_to_use = list(range(3, 4)) # FIXME
    # boxes_to_use.reverse()

    Ds[file] = {}
    D_uncs[file] = {}

    data = common.load(f'box_counting/data/counted_{file}.npz')
    # data = common.load(f'data/counted_driftremoved_{phi}.npz')
    N2_mean = data['N2_mean']
    N2_std  = data['N2_std']
    # N_stats = data['N_stats']
    phi     = data['pack_frac']
    sigma   = data['particle_diameter']
    time_step    = data['time_step']
    depth_of_field = data.get('depth_of_field')

    box_sizes = data['box_sizes']
    sep_sizes = data['sep_sizes']
    # N_mean    = N_stats[:, 1]
    # N_var     = N_stats[:, 2]
    # num_boxes_used = N_stats[:, 5]

    num_timesteps = N2_mean.shape[1]
    num_boxes     = N2_mean.shape[0]
    t_all = np.arange(0, num_timesteps) * time_step

    t_theory = np.logspace(np.log10(t_all[1] / 5), np.log10(t_all[-1]))
    # print(t_theory,np.log10(t[1] / 10), np.log10(t.max()))


    # D0 = { # countoscope paper, table 1
    #     0.02: 0.0416,
    #     0.34: 0.0310,
    #     0.66: 0.0175
    # }[phi]

    # if mode == 'strips':
    #     D0 = D0 / 2

    # for box_size_index in boxes_to_use:
    for box_size_index in range(len(box_sizes)):
    # for box_size_index in [1]:
    # for L in [2**e for e in range(-2, 7)]:
        L = box_sizes[box_size_index]
        t = np.copy(t_all)

        delta_N_sq = N2_mean[box_size_index, :]
        # t = np.arange(0, len(delta_N_sq))[1:]/2
        # delta_N_sq = delta_N_sq # [1:] is because the point at t=0 msd=0 plots weirdly
        
        # D0 = 0.038 # Bare diffusion coefficient in um^2/s -- short time?
        # # D0 = D0 * 2.2

        # if phi == 0.66:   
        #     D0 = D0 / 2.2 / 1.2
        #     pass

        # first_grad = ]
        f = lambda tau: np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2
        
        L_2 = L
        
        # N2_func = lambda t, D0: 8/np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t / L**2) # countoscope eq. 3
        # N2_func_full = lambda t, D0: 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2) * f(4*D0*t/L_2**2)) # countoscope eq. 2, countoscope overleaf doc

        # fit_func = N2_func_full
        # popt, pcov = scipy.optimize.curve_fit(fit_func, t[0:LOWTIME_FIT_END], N2_mean[box_size_index, 0:LOWTIME_FIT_END])
        # D0 = popt[0]
        # r2 = common.r_squared(N2_mean[box_size_index, 0:LOWTIME_FIT_END], fit_func(t[0:LOWTIME_FIT_END], D0))
        # fit to whole thing
        if depth_of_field:
            N2_theory = lambda t, D, N: common.N2_nointer_3D(t, D, N, L, L, depth_of_field)
        else:
            N2_theory = lambda t, D, N: common.N2_nointer_2D(t, D, N, L, L)
        fitting_points = np.unique(np.round(10**np.linspace(0, np.log10(t.max()/2))).astype('int'))
        print(fitting_points)
        popt, pcov = scipy.optimize.curve_fit(N2_theory, t[fitting_points], delta_N_sq[fitting_points], maxfev=2000)
        
        Ds[file][L] = popt[0]
        D_uncs[file][L] = np.sqrt(pcov[0][0])

        if collapse_x:
            t /= L**2
        if collapse_y:
            delta_N_sq /= N_var[box_size_index]

        # t /= t[-1]

        #, r^2={r2:.2f}
        # D_str += f'±{np.sqrt(pcov[0][0]):.3f}'

        # ax.plot(t_theory, N2_func_full(t_theory, D0), color='black', zorder=5, linestyle='dotted', linewidth=1, label='sFDT (no inter.)' if box_size_index==0 else None)

        # ax.hlines(2*N_mean[box_size_index], t.min(), t.max(), color='grey', linewidth=1, label=r'$2 \langle N \rangle$' if box_size_index==0 else None)
        # ax.hlines(2*N_var [box_size_index], t.min(), t.max(), linestyles='dashed', color='grey', linewidth=1, label=r'$\mathrm{Var}(N)$' if box_size_index==0 else None)


        # N2_theory = 2 * N_mean[box_size_index] * (1 - f(4*D0*t/L**2)**2) # countoscope eq. 2
        # N2_theory_lowtime = 4 / np.sqrt(np.pi) * N_mean[box_size_index] * np.sqrt(D0 * t_theory) * (L + L_2) / (L * L_2)
        # ax.plot(t_theory[:LOWTIME_FIT_END], N2_theory_lowtime[:LOWTIME_FIT_END], linestyle='dashed', linewidth=1, color='black', label='sFDT (no inter.) low time' if box_size_index==0 else None)

        # p1, p2 = plateaus.calc_plateaus_for_L(sigma, phi, L)
        # ax.hlines(p1, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1, label='plateaus')
        # ax.hlines(p2, t.min(), t.max(), linestyles='dashed', color=exp_plot[0].get_color(), linewidth=1)

        # N2_theory_interactions = 2 * N_var[box_size_index] * sDFT_interactions.sDFT_interactions(L, t_theory, phi, D0, sigma)# * 10
        # ax.plot(t_theory, N2_theory_interactions, color='black', linestyle='dotted', linewidth=1, zorder=3, label='sFDT (w/ inter.)' if box_size_index==0 else None)

        label = rf'$L = {L}\mathrm{{\mu m}}$'
        label += f', {file}'
        # label += f', D={D0:.3f}'
        # label += f', $sep = {sep_sizes[box_size_index]:.1f}\mathrm{{\mu m}}$'
        # label += f', $n = {num_boxes_used[box_size_index]:.0f}$'
        
        if frac := data.get('data_fraction'):
            label += f', {frac:.3f} used'

        if collapse_x or collapse_y:
            markersize = 2
        else:
            markersize = 2
        ax.plot(t[1:], delta_N_sq[1:], label=label, color=color, linestyle='none', marker='o', markersize=markersize, markeredgecolor='none')
        # ax.plot(t_theory[1:], N2_theory(t_theory, *popt)[1:], color='black', linewidth=1, label='sDFT (no inter.)' if box_size_index==0 else None)
        

    titles.append(f'{file}, $\phi_\mathrm{{calc}}={phi:.3f}$, $\sigma={sigma}$')

    # ax.legend(loc='lower right', fontsize=8)
legend = ax.legend(fontsize=6, loc='lower right')
for handle in legend.legend_handles:
    handle.set_markersize(6.0)
ax.semilogy()
ax.semilogx()
xlabel = '$t/L^2$' if collapse_x else '$t$'
ylabel = r'$\Delta N^2(t) / \rangle N \langle$' if collapse_y else '$\Delta N^2(t)$'
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(', '.join(titles))
ax.set_title(titles[0])

# if collapse_x:
#     ax.set_ylim(1e-1, 4e0)
# else:
#     ax.set_ylim(1e-4, 2e-2)
# print('NAUGHTY YLIM')

fig.tight_layout()
names = '_'.join(sys.argv[1:])
common.save_fig(fig, f'box_counting/figures_png/msd_combined_{names}.png', dpi=300)


sigma_calced = data.get('particle_diameter_calced')

temp = []

# D graph
fig, ax = plt.subplots(1, 1)
for file_i, file in enumerate(files):
    # print(Ds[file])
    for L, D in Ds[file].items():
        color = matplotlib.cm.afmhot(0.3+np.log(L)/6)
        # print(0.3+np.log(L)/6)
        # ax.scatter([file], [D], color=color, label=fr'$L={L:.2f}\mathrm{{\mu m}}$' if file_i==0 else None)
        D_unc = D_uncs[file][L]
        
        if D_unc * 1 > D:
            continue

        temp.append(D*sigma_calced)

        ax.errorbar([file], [D], yerr=D_unc, marker='o', color=color, label=fr'$L={L:.2f}\mathrm{{\mu m}}$')
# ax.set_ylim(0, 2)
# ax.semilogy()
plt.xticks(rotation=45)
# plt.xticks(rotation=45, ha="right")
ax.set_ylabel('$D$')

print(np.mean(temp), np.std(temp))

# remove duplicate legends and sort
# https://stackoverflow.com/a/13589144/1103752
# https://stackoverflow.com/a/46160465/1103752
handles, labels = ax.get_legend_handles_labels()
zipped_list = sorted(list(zip(labels, handles)), key=lambda t: t[0])
by_label = dict(zipped_list) # removes duplicates
ax.legend(by_label.values(), by_label.keys(),
          fontsize=8,
          bbox_to_anchor=(1, 1), loc='upper left')
# ax.legend(bbox_to_anchor=(1, 1))

common.save_fig(fig, f'box_counting/figures_png/D_combined_{names}.png')