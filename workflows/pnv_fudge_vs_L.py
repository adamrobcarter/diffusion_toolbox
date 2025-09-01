import common
import particle_detection.add_drift_periodic
import particle_linking.add_drift_periodic
import box_counting.calc_pnv
import box_counting.show_pnv

import matplotlib.pyplot as plt
import numpy as np
import particle_linking.link

# vs = [0.25, 0.5, 1]
# vs = [0.5, 1]
# vs = [1]

v_slice = 0
phi_slice = 0

# file_base = 'sim_nointer'
# L = '640'
# phis = ['001', '010', '030', '050']

file_base = 'sim_nohydro'
L = '320_short'
# phis = ['001', '010', '020', '030', '040', '050']
# phis = ['001', '010', '020', '030', '040']
phis = [0.005, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25]
# phis = ['001', '050']

# N0_sources = ['mean', 'var', 'special']
# N0_sources = ['N0S0']
N0_sources = ['var', 'N0S0']

# v_i = 0
v_true = 1

aspects = [0.2, 1, 5]

F_MULT_NONE = 0
F_MULT_VAR = 1
F_MULT = F_MULT_VAR

## do precaculation
for phi_i, phi in enumerate(phis):
    # particle_linking.link.go(f'{file_base}_{phi}_L{L}')

    # for v_i, v in enumerate(vs):
    # particle_detection.add_drift_periodic.go(f'{file_base}_{phi}_L{L}', 'const', v_true)
    
    # linked simple method
    # particle_linking.add_drift_periodic.go(f'{file_base}_{phi}_L{L}', 'const', v)
    # data = common.load(f'particle_linking/data/trajs_{file_base}_{phi}_L{L}_drifted_const_v{v}.npz')
    # drift = common.find_drift(data['particles'], data.get('dimension', 2))
    # print('drift', drift)

    # pnv method
    # for aspect in aspects:
    #     box_counting.calc_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v_true}', frame_deltas=range(0, 25), aspect=aspect)
    pass

# for N0_source in N0_sources:
for aspect in aspects:
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    fudge_mean = np.full((len(phis)), np.nan)
    fudge_std  = np.full((len(phis)), np.nan)

    pnv_vs = np.full((len(phis), 20), np.nan)
    pnv_Ls = None

    for phi_i, phi in enumerate(phis):
        # for v_i, v in enumerate(vs):

        # linked simple method
        # data = common.load(f'particle_linking/data/trajs_{file_base}_{phi}_L{L}_drifted_const_v{v}.npz')
        # drift = common.find_drift(data['particles'], data.get('dimension', 2))
        
        # pnv method
        file = f'{file_base}_{phi}_L{L}_drifted_const_v{v_true}_aspect{aspect}'
        data = common.load(f'box_counting/data/pnv_{file}.npz')
        time_step = data['time_step']
        N1N2 = data['N1N2'] # dimensions are box size, t
        # N1N2_std = data['N1N2_std']
        # N1N2_std = data['N1N2_per_box_unc'].mean(axis=(1, 2)) # average over all boxes
        N1N2_err = np.zeros_like(N1N2)
        L1 = data['box_sizes_x']
        L2 = data['box_sizes_y']
        N0 = data['N_mean']
        N_var = data['N_var']
        N1N2mN1N2 = data['N1N2mN1N2']
        particle_diameter = data['particle_diameter']

        F_N0 = N1N2[:, 1] * L1 / (time_step * v_true)
        F_N0_unc = N1N2_err[:, 1] * L1 / (time_step * v_true)

        if F_MULT == F_MULT_NONE:
            F = F_N0 / N0
            F_unc = F_N0_unc / N0
        elif F_MULT == F_MULT_VAR:
            F = F_N0 / N_var
            F_unc = F_N0_unc / N_var


        pnv_vs[phi_i, :len(F)] = F
        pnv_Ls = L1

        # pnv_V, pnv_L, time_step = box_counting.show_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v_true}', N0_source=N0_source)
        # pnv_vs[phi_i, :len(pnv_V)] = pnv_V
        # if phi_i > 0:
        #     assert np.all(pnv_Ls == pnv_Ls)
        # pnv_Ls = pnv_L

        # print()

        # print('min L', v_true * time_step)
        # for L_i in range(len(pnv_Ls)):
        #     if pnv_Ls[L_i] < v_true * time_step:
        #         print('skipped')
        #         continue
        to_remove = L1 < v_true * time_step
        # print(L1, L2)
        print(f'length threshold = v dt = {v_true * time_step:.2g}um (num: {to_remove.sum()})')
        L1    = L1   [~to_remove]
        L2    = L2   [~to_remove]
        F     = F    [~to_remove]
        F_unc = F_unc[~to_remove]

        ax.errorbar(np.sqrt(L1*L2)/particle_diameter, F, yerr=F_unc, marker='o', label=f'$\phi={phi}$', color=common.colormap_cool(phi_i, 0, len(phis)))
        # ax_phi.plot(np.sqrt(L1*L2), F, marker='o', label=f'$L_x={pnv_Ls[L_i]/particle_diameter:.2g}\sigma$', color=common.colormap(L_i, 0, len(pnv_Ls)))
# 
    # ax_phi.errorbar(phis_value, fudge_mean, yerr=fudge_std, label='PNV', marker='o', linestyle='none')

    ax.set_xlabel(f'$\sqrt{{L_xL_y}}/\sigma$')
    if F_MULT == F_MULT_NONE:
        ax.set_ylabel('$F(\phi, L, \sigma)$')
    else:
        ax.set_ylabel('$F(\phi, L, \sigma) N_0 / \mathrm{{Var}}(N)$')
    # ax_phi.set_ylim(0, 2)

    # if file_base == 'sim_nointer':
    #     ax_phi.set_ylim(0, 1.5*v_true)

    ax.legend(fontsize=7)

    ax.grid(alpha=0.5)

    ax.set_ylim(0, 1.1)
    # ax.semilogx()

    fig.suptitle({'sim_nointer': 'no interactions', 'sim_nohydro': 'steric interactions'}[file_base] + f', $L_x/L_y={aspect}$')

    filename = f'fudge_factors_vs_L_{file_base}_aspect{aspect}'
    if F_MULT == F_MULT_VAR:
        filename += '_var'
    common.save_fig(fig, f'box_counting/figures_png/{filename}.png')