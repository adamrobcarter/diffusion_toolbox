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


VAR = '\mathrm{{Var}}'

for file in common.files_from_argv('box_counting/data/', 'pnv_'):
    data = common.load(f'box_counting/data/pnv_{file}_special.npz')
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

    to_remove = L1[1:-3:3] < v_true * time_step
    print(f'length threshold = v dt = {v_true * time_step:.2g}um (num: {to_remove.sum()})')

    x = np.sqrt(L1*L2)/particle_diameter

    print('plotting first 3')
    ax.plot(x[1:-3:3][~to_remove], F_N0[1:-3:3][~to_remove], marker='o', label='PNV')
    ax.plot(x[1:-3:3], N0[1:-3:3],   marker='o', label='$N_0$')
    ax.plot(x[1:-3:3], N_var[1:-3:3],   marker='o', label=fr'$\sqrt{{{VAR}(N_1){VAR}(N_2)}}$')
    
    # N1N2mN1N2 = - np.nanmean(N1, axis=(1, 2, 3)) * np.nanmean(N2, axis=(1, 2, 3)) - np.nanmean(N1 * N2, axis=(1, 2, 3))
    # print('N1N2mN1N2', N1N2mN1N2)
    # ax.plot(x, N1N2mN1N2,   marker='o', label=r'$\langle N_1 \rangle \langle N_2 \rangle - \langle N_1 N_2 \rangle$')
        # ax_phi.plot(np.sqrt(L1*L2), F, marker='o', label=f'$L_x={pnv_Ls[L_i]/particle_diameter:.2g}\sigma$', color=common.colormap(L_i, 0, len(pnv_Ls)))
# 
    # ax_phi.errorbar(phis_value, fudge_mean, yerr=fudge_std, label='PNV', marker='o', linestyle='none')

    print('assertions')
    assert np.any(N1 != 0)
    assert np.any(N2 != 0)

    # sophie2 = (np.nanvar(N1[1::3, :, :, :], axis=(1, 2, 3))+np.nanvar(N2[1::3, :, :, :], axis=(1, 2, 3)))/2 - np.nanvar(N1[1::3, :, :, :]*N2[1::3, :, :, :], axis=(1, 2, 3))
    # ax.plot(x[1::3], sophie2, marker='o', label=rf'$({VAR}(N_1)+{VAR}(N_2))/2 - {VAR}(N_1N_2)$')
    
    # sophie3 = np.sqrt(np.nanvar(N1[1::3, :, :, :], axis=(1, 2, 3))*np.nanvar(N2[1::3, :, :, :], axis=(1, 2, 3))) - np.nanvar(N1[1::3, :, :, :]*N2[1::3, :, :, :], axis=(1, 2, 3))
    # print(sophie3)
    # ax.plot(x[1::3], sophie3, marker='o', label=rf'$\sqrt{{{VAR}(N_1){VAR}(N_2)}} - {VAR}(N_1N_2)$')
    
    print('p')
    Var12 = np.nanmean(N1, axis=(1, 2, 3))*np.nanmean(N2, axis=(1, 2, 3)) - np.nanmean(N1*N2, axis=(1, 2, 3))
    print('p')
    Var12_double = Var12[1::3] # only the real boxes
    print('p')
    Var12_double = Var12_double[3:] # shifts so we now have the one of 2L
    print('p')
    Var12_double = np.concatenate([Var12_double, [np.nan, np.nan, np.nan]])
    print('p')

    d_Var12 = Var12[2::3] - Var12[0::3]
    d_L = L1[2::3] - L1[0::3]

    print('Var12', Var12)
    print('d_Var12', d_Var12)
    print('d_L', d_L)
    print('rho var', data['density'] * Var12_double )
    print('d/d', d_Var12/d_L)
    # print('theyare', data['density'], Var12_double, d_Var12, d_L)

    # alpha = data['density'] * Var12_double + d_Var12/d_L
    alpha = Var12_double + d_Var12/d_L
    alpha2 = data['density'] * Var12_double + d_Var12/d_L

    ax.plot(L1[1::3]/particle_diameter, Var12_double, marker='o', label=fr'$-{VAR}_{{12}}(2L_1,2L_1)$')
    ax.plot(L1[1::3]/particle_diameter, d_Var12/d_L, marker='o', label=fr'$-\partial{VAR}_{{12}}(L_1,L_1)/\partial L_1$')
    ax.plot(L1[1::3]/particle_diameter, alpha, marker='o', label=fr'$-{VAR}_{{12}}(2L_1,2L_1)-\partial{VAR}_{{12}}(L_1,L_1)/\partial L_1$')

    ax.set_xlabel(f'$\sqrt{{L_xL_y}}/\sigma$')
    ax.set_ylabel('$F(\phi, L, \sigma)N_0$')
    # ax_phi.set_ylim(0, 2)

    # if file_base == 'sim_nointer':
    #     ax_phi.set_ylim(0, 1.5*v_true)

    ax.legend(fontsize=7)

    ax.grid(alpha=0.5)

    # ax.set_ylim(0, 1.1)
    ax.semilogy()
    ax.semilogx()

    # fig.suptitle(file)
    ax.set_title(f'$\phi={phi:.3f}$, $L_x/L_y={np.mean(L1/L2)}$')

    common.save_fig(fig, f'box_counting/figures_png/fudge_factors_vs_L_special_{file}.png')