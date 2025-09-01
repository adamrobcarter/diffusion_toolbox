import common
import particle_detection.add_drift_periodic
import particle_linking.add_drift_periodic
import box_counting.calc_pnv
import box_counting.show_pnv

import matplotlib.pyplot as plt
import numpy as np
import particle_linking.link

vs = [0.25, 0.5, 1]
# vs = [0.5, 1]
# vs = [1]

v_slice = 0
phi_slice = 0

# file_base = 'sim_nointer'
# L = '640'
# phis = ['001', '010', '030', '050']

file_base = 'sim_nohydro'
L = '320_short'
phis = ['001', '010', '020', '030', '040', '050']
# phis = ['001', '050']

# N0_sources = ['mean', 'var', 'special']
# N0_sources = ['N0S0']
N0_sources = ['var', 'N0S0']

v_i = 0
v_true = vs[v_i]

aspects = [0.2, 1, 5]

## do precaculation
for phi_i, phi in enumerate(phis):
    # particle_linking.link.go(f'{file_base}_{phi}_L{L}')

    # for v_i, v in enumerate(vs):
    # particle_detection.add_drift_periodic.go(f'{file_base}_{phi}_L{L}', 'const', v)
    
    # linked simple method
    # particle_linking.add_drift_periodic.go(f'{file_base}_{phi}_L{L}', 'const', v)
    # data = common.load(f'particle_linking/data/trajs_{file_base}_{phi}_L{L}_drifted_const_v{v}.npz')
    # drift = common.find_drift(data['particles'], data.get('dimension', 2))
    # print('drift', drift)

    # pnv method
    for aspect in aspects:
        box_counting.calc_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v_true}', frame_deltas=[0, 1], aspect=aspect)
    pass

# for N0_source in N0_sources:
for aspect in aspects:
    fig, ax_phi = plt.subplots(1, 1, figsize=(4, 4))

    fudge_mean = np.full((len(phis)), np.nan)
    fudge_std  = np.full((len(phis)), np.nan)
    phis_value = np.array([int(phi)/100 for phi in phis])

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
        L1 = data['box_sizes_x']
        L2 = data['box_sizes_y']
        N0 = data['N_mean']
        N1N2mN1N2 = data['N1N2mN1N2']
        print('N0', N0)
        print('N1N2mN1N2', N1N2mN1N2)
        particle_diameter = data['particle_diameter']

        N1N2_frame1 = N1N2[:, 1]

        F = N1N2_frame1 / N0 * L1 / (time_step * v_true)

        pnv_vs[phi_i, :len(F)] = F
        pnv_Ls = L1

        # pnv_V, pnv_L, time_step = box_counting.show_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v_true}', N0_source=N0_source)
        # pnv_vs[phi_i, :len(pnv_V)] = pnv_V
        # if phi_i > 0:
        #     assert np.all(pnv_Ls == pnv_Ls)
        # pnv_Ls = pnv_L

        # print()

    print('min L', v_true * time_step)
    for L_i in range(len(pnv_Ls)):
        if pnv_Ls[L_i] < v_true * time_step:
            print('skipped')
            continue



        ax_phi.plot(phis_value, pnv_vs[:, L_i], marker='o', label=f'$L_x={pnv_Ls[L_i]/particle_diameter:.2g}\sigma$', color=common.colormap(L_i, 0, len(pnv_Ls)))
        ax_phi.plot(phis_value, pnv_vs[:, L_i], marker='o', label=f'$L_x={pnv_Ls[L_i]/particle_diameter:.2g}\sigma$', color=common.colormap(L_i, 0, len(pnv_Ls)))

    # ax_phi.errorbar(phis_value, fudge_mean, yerr=fudge_std, label='PNV', marker='o', linestyle='none')

    ax_phi.set_xlabel(f'packing fraction')
    ax_phi.set_ylabel('$F(\phi, L, \sigma)$')
    # ax_phi.set_ylim(0, 2)

    # if file_base == 'sim_nointer':
    #     ax_phi.set_ylim(0, 1.5*v_true)

    ax_phi.legend(fontsize=7)

    ax_phi.grid(alpha=0.5)

    ax_phi.set_ylim(0, 1.1)

    fig.suptitle({'sim_nointer': 'no interactions', 'sim_nohydro': 'steric interactions'}[file_base] + f', $L_x/L_y={aspect}$')

    common.save_fig(fig, f'box_counting/figures_png/fudge_factors_new_{file_base}_aspect{aspect}.png')