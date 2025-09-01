import common
import matplotlib.pyplot as plt
import numpy as np
import scipy

fig, ax = plt.subplots(1, 1)

phis = []
drifts = []
drift_uncs = []
# zetas = []

# for file in common.files_from_argv('particle_linking/data', 'trajs_'):
#     data = common.load(f'particle_linking/data/trajs_{file}.npz')
#     particle_diameter = data['particle_diameter']   
#     drift_xyz, drift_xyz_unc = common.find_drift(data['particles'], data.get('dimension', 2))
#     # z_avg = data['particles'][:, 2].mean()
#     # z_zeta = z_avg / (particle_diameter/2) - 1

#     # print('drift', drift_xyz)
#     phis.append(data['pack_frac'])
#     drifts.append(drift_xyz[0]) # um/s
#     drift_uncs.append(drift_xyz_unc[0]) # um/s

#     # theta = data['theta']
#     # zetas.append(z_zeta)

# ax.scatter(phis, drifts/np.sin(10*np.pi/180), label='sim')

# theta_theory = np.linspace(min(phis), max(phis))

# zetas = np.array(zetas)
# phis = np.array(phis)

eleanor_phis = np.array([0.068, 0.216, 0.256, 0.299, 0.307, 0.357, 0.453, 0.545, 0.545])
eleanor_vs = np.array([0.692, 0.875, 0.862, 0.898, 0.94, 0.982, 1.12, 1.17, 1.19])

ax.scatter(eleanor_phis, eleanor_vs, label='exp')

ax2 = ax.twinx()
Sk0s = [common.S_k_zero(phi) for phi in eleanor_phis]
ax2.scatter(eleanor_phis, Sk0s)

ax.legend()
ax.set_xlabel(r'$\phi$')
ax.set_ylabel(r'$v/\mathrm{sin}\theta$')

common.save_fig(fig, f'particle_linking/figures_png/drift_vs_phi_{file}.png')