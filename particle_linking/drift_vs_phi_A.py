import common
import matplotlib.pyplot as plt
import numpy as np
import scipy
import visualisation.Ds_overlapped
import re

fig, ax = plt.subplots(1, 1)

phis = []
As = []
drift_uncs = []
# zetas = []

# python -m particle_linking.drift_vs_phi sim_hydro_*_L640_theta10 sim_nohydro_*_L640_theta10_unwrap

if False:
    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particle_diameter = data['particle_diameter']   
        pack_frac_given = data['pack_frac_given']
        drift_xyz, drift_xyz_unc = common.find_drift(data['particles'], data.get('dimension', 2))
        # z_avg = data['particles'][:, 2].mean()
        # z_zeta = z_avg / (particle_diameter/2) - 1

        # need to get the theta0 MSD. this is more complicated because we changed naming convention in the middle
        # also should it really be the theta0 D or the theta10 D?
        hydropart = 'nohydro' if 'nohydro' in file else 'hydro'
        # D0_file = f'sim_{hydropart}_{pack_frac_given}_L640_theta0_1h_unwrap'
        D0_file = file.replace('theta10', 'theta0')
        D0_file = re.sub('hydro_0\.[0-9]+', 'hydro_0.02', D0_file)
        # sim_hydro_0.08_L640_theta0_1h_unwrap
        D, _, _ = visualisation.Ds_overlapped.get_D0(D0_file)

        # print('drift', drift_xyz)
        phis.append(data['pack_frac'])
        delta_M_g = 0.0592 # assumed, taken from inputfile. once you switch to multiblobs you can do this properly
        kBT = 0.0041419464 # taken from inputfile. once you switch to multiblobs you can do this properly
        norm = np.sin(10*np.pi/180) * D * delta_M_g / kBT
        norm /= 2 # DIRTY HACK
        A = drift_xyz[0] / norm
        A_unc = drift_xyz_unc[0] / norm
        As.append(A) # um/s
        drift_uncs.append(A_unc) # um/s

        # theta = data['theta']
        # zetas.append(z_zeta)

    phis = np.array(phis)
    As = np.array(As)

    ax.scatter(phis, As, label='sim')
    print(phis)
    print(As)

# theta_theory = np.linspace(min(phis), max(phis))

# zetas = np.array(zetas)
# phis = np.array(phis)

# eleanor_phis = np.array([0.068, 0.216, 0.256, 0.299, 0.307, 0.357, 0.453, 0.545, 0.545])
# eleanor_vs = np.array([0.692, 0.875, 0.862, 0.898, 0.94, 0.982, 1.12, 1.17, 1.19])
eleanor_phis = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1])
eleanor_As = np.array([0.92, 0.89, 0.95, 1.01, 1.06, 1.11, 1.18])

ax.scatter(eleanor_phis, eleanor_As, label='exp')


ax2 = ax.twinx()
Sk0s = [common.S_k_zero(phi) for phi in eleanor_phis]
ax2.plot(eleanor_phis, Sk0s, color='tab:red', )

ax.legend()
ax.set_xlabel(r'$\phi$')
ax.set_ylabel(r'$A = v/(\mathrm{{sin}}(\theta) D(\phi \rightarrow 0) \Delta m g / k_B T)$')

file = ''
common.save_fig(fig, f'particle_linking/figures_png/drift_vs_phi_A_{file}.png')