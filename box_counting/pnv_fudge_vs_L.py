import common
import particle_detection.add_drift_periodic
import particle_linking.add_drift_periodic
import box_counting.calc_pnv
import box_counting.show_pnv

import matplotlib.pyplot as plt
import numpy as np
import particle_linking.link


# N0_sources = ['mean', 'var', 'special']
# N0_sources = ['N0S0']
N0_sources = ['var', 'N0S0']

for file in common.files_from_argv('box_counting/data/', 'pnv_'):
    data = common.load(f'box_counting/data/pnv_{file}.npz')
    phi = data['pack_frac']
    v_true = data['velocity_multiplier']
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    time_step = data['time_step']
    N1N2 = data['N1N2'] # dimensions are box size, t
    # N1N2_std = data['N1N2_std']
    # N1N2_std = data['N1N2_per_box_unc'].mean(axis=(1, 2)) # average over all boxes
    N1N2_err = np.zeros_like(N1N2)
    L1        =  data['box_sizes_x']
    L2        =  data['box_sizes_y']
    N0        =  data['N_mean']
    N_var     =  data['N_var']
    particle_diameter = data['particle_diameter']
    N1 = data['counts_N1']
    N2 = data['counts_N2']

    F_N0 = N1N2[:, 1] * L1 / (time_step * v_true)

    to_remove = L1 < v_true * time_step
    print(f'length threshold = v dt = {v_true * time_step:.2g}um (num: {to_remove.sum()})')

    x = np.sqrt(L1*L2)/particle_diameter
    ax.plot(x[~to_remove], F_N0[~to_remove], marker='o', label='PNV')
    ax.plot(x, N0,   marker='o', label='$N_0$')
    ax.plot(x, N_var,   marker='o', label=r'$\sqrt{{{var}(N_1){var}(N_2)}}$')
    
    # N1N2mN1N2 = - np.nanmean(N1, axis=(1, 2, 3)) * np.nanmean(N2, axis=(1, 2, 3)) - np.nanmean(N1 * N2, axis=(1, 2, 3))
    # print('N1N2mN1N2', N1N2mN1N2)
    # ax.plot(x, N1N2mN1N2,   marker='o', label=r'$\langle N_1 \rangle \langle N_2 \rangle - \langle N_1 N_2 \rangle$')
        # ax_phi.plot(np.sqrt(L1*L2), F, marker='o', label=f'$L_x={pnv_Ls[L_i]/particle_diameter:.2g}\sigma$', color=common.colormap(L_i, 0, len(pnv_Ls)))
# 
    # ax_phi.errorbar(phis_value, fudge_mean, yerr=fudge_std, label='PNV', marker='o', linestyle='none')

    var = r'\mathrm{{Var}}'
    # assert np.any(N1 != 0)
    # assert np.any(N2 != 0)
    # sophie2 = (np.nanvar(N1, axis=(1, 2, 3))+np.nanvar(N2, axis=(1, 2, 3)))/2 - np.nanvar(N1*N2, axis=(1, 2, 3))
    # ax.plot(x, sophie2, marker='o', label=rf'$({var}(N_1)+{var}(N_2))/2 - {var}(N_1N_2)$')
    
    # sophie3 = np.sqrt(np.nanvar(N1, axis=(1, 2, 3))*np.nanvar(N2, axis=(1, 2, 3))) - np.nanvar(N1*N2, axis=(1, 2, 3))
    # print(sophie3)
    # ax.plot(x, sophie3, marker='o', label=rf'$\sqrt{{{var}(N_1){var}(N_2)}} - {var}(N_1N_2)$')
    
    ax.set_xlabel(fr'$\sqrt{{L_xL_y}}/\sigma$')
    ax.set_ylabel(r'$F(\phi, L, \sigma)N_0$')
    # ax_phi.set_ylim(0, 2)

    # if file_base == 'sim_nointer':
    #     ax_phi.set_ylim(0, 1.5*v_true)

    ax.legend(fontsize=7)

    ax.grid(alpha=0.5)

    # ax.set_ylim(0, 1.1)
    ax.semilogy()
    ax.semilogx()

    # fig.suptitle(file)
    ax.set_title(fr'$\phi={phi:.3f}$, $L_x/L_y={np.mean(L1/L2)}$')

    common.save_fig(fig, f'box_counting/figures_png/fudge_factors_vs_L_{file}.png')