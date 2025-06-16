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
# phis = ['001', '010', '020', '030', '040', '050']
phis = ['020', '040']

# N0_sources = ['mean', 'var', 'special']
N0_sources = ['N0S0']

for N0_source in N0_sources:
    fig, (ax_v, ax_phi) = plt.subplots(1, 2, figsize=(7, 4))

    fudge_mean = np.full((len(vs), len(phis)), np.nan)
    fudge_std  = np.full((len(vs), len(phis)), np.nan)
    phis_value = np.array([int(phi)/100 for phi in phis])

    for phi_i, phi in enumerate(phis):
        # particle_linking.link.go(f'{file_base}_{phi}_L{L}')

        for v_i, v in enumerate(vs):
            # linked simple method
            # print('v', v)
            # particle_linking.add_drift_periodic.go(f'{file_base}_{phi}_L{L}', 'const', v)
            # data = common.load(f'particle_linking/data/trajs_{file_base}_{phi}_L{L}_drifted_const_v{v}.npz')
            # drift = common.find_drift(data['particles'], data.get('dimension', 2))
            # print('drift', drift)

            # pnv method
            # particle_detection.add_drift_periodic.go(f'{file_base}_{phi}_L{L}', 'const', v)
            # box_counting.calc_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v}')
            mean, std = box_counting.show_pnv.go(f'{file_base}_{phi}_L{L}_drifted_const_v{v}', N0_source=N0_source)

            # assert 0.8 < mean / v < 1.2
            print()

            fudge_mean[v_i, phi_i] = mean / v
            fudge_std [v_i, phi_i] = std  / v

    ax_v.errorbar(vs, fudge_mean[:, phi_slice], yerr=fudge_std[:, phi_slice], label='phi=0.10', marker='o', linestyle='none')
    ax_phi.errorbar(phis_value, fudge_mean[v_slice, :], yerr=fudge_std[v_slice, :], label='v=0.50', marker='o', linestyle='none')

    ax_v  .set_xlabel(fr'added velocity ($\mathrm{{\mu m s^{{-1}}}}$) ($\phi={float(phis[phi_slice])/100}$)')
    ax_phi.set_xlabel(f'packing fraction ($v={vs[v_slice]}$)')
    ax_v  .set_ylabel('$v_{PNV} / v_{true}$')
    ax_phi.set_ylabel('$v_{PNV} / v_{true}$')
    ax_phi.set_ylim(0, 2)

    fig.suptitle({'sim_nointer': 'no interactions', 'sim_nohydro': 'steric interactions'}[file_base] + ', ' + 'N0: ' + N0_source)

    common.save_fig(fig, f'box_counting/figures_png/fudge_factors_{file_base}_{N0_source}.png')