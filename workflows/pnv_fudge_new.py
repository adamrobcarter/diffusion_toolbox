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
v = vs[v_i]

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
    # box_counting.calc_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v}', frame_deltas=range(0, 25))
    pass

for N0_source in N0_sources:
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
        pnv_V, pnv_L, dt = box_counting.show_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v}', N0_source=N0_source)
        pnv_vs[phi_i, :len(pnv_V)] = pnv_V
        if phi_i > 0:
            assert np.all(pnv_Ls == pnv_Ls)
        pnv_Ls = pnv_L

        print()

    print('min L', v * dt)
    for L_i in range(len(pnv_Ls)):
        if pnv_Ls[L_i] < v * dt:
            print('skipped')
            continue

        ax_phi.scatter(phis_value, pnv_vs[:, L_i], label=f'PNV measured L={pnv_Ls[L_i]:.1f}', color=common.colormap(L_i, 0, len(pnv_Ls)))

    # ax_phi.errorbar(phis_value, fudge_mean, yerr=fudge_std, label='PNV', marker='o', linestyle='none')

    ax_phi.set_xlabel(f'packing fraction')
    ax_phi.set_ylabel('$v$ ($\mathrm{\mu m s^{-1}}$)')
    # ax_phi.set_ylim(0, 2)

    if file_base == 'sim_nointer':
        ax_phi.set_ylim(0, 1.5*v)

    ax_phi.hlines(v, *ax_phi.get_xlim(), label='true velocity')
    ax_phi.legend(fontsize=7)

    ax_phi.grid(alpha=0.5)

    fig.suptitle({'sim_nointer': 'no interactions', 'sim_nohydro': 'steric interactions'}[file_base] + ', ' + 'N0: ' + N0_source)

    common.save_fig(fig, f'box_counting/figures_png/fudge_factors_new_{file_base}_{N0_source}.png')