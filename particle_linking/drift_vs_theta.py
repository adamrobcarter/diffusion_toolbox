import common
import matplotlib.pyplot as plt
import numpy as np
import scipy

def go(files, ax):
    thetas = []
    drifts = []
    drift_uncs = []
    zetas = []

    for file in files:
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particle_diameter = data['particle_diameter']   
        drift_xyz, drift_xyz_unc = common.find_drift(data['particles'], data.get('dimension', 2))
        z_avg = data['particles'][:, 2].mean()
        z_zeta = z_avg / (particle_diameter/2) - 1

        # print('drift', drift_xyz)
        thetas.append(data['theta'])
        drifts.append(-drift_xyz[0]) # um/s
        drift_uncs.append(drift_xyz_unc[0]) # um/s
        zetas.append(z_zeta)
        n_blobs = data['num_blobs']

    ax.errorbar(thetas, drifts, yerr=drift_uncs, linestyle='none', marker='o', label=fr'simulation ($n_\mathrm{{blobs}}={n_blobs}$)')
    for i in range(len(thetas)):
        print(f'θ={thetas[i]}, v={drifts[i]:.5f}um/s')


    theta_theory = np.linspace(min(thetas), max(thetas))

    g = scipy.constants.g # m/s^2
    g *= 1e6 # um/s^2
    delta_rho = data['delta_rho'] # kg/um^3
    eta = data['eta'] # Pa s == kg/m/s
    eta *= 1e-6 # convert to kg/um/s
    # from https://www.mdpi.com/2076-3263/10/2/69
    theory_niiya = - 2/9 * g * delta_rho / eta * np.sin(theta_theory * np.pi/180) * (particle_diameter/2)**2
    # units = um/s^2 * kg/um^3 / kg/um/s * um^2 = um/s
    # ax.plot(theta_theory, theory_niiya, label='Niiya 2019')


    # the following from https://www.sciencedirect.com/science/article/pii/S0301932202000770
    zetas = np.array(zetas)
    thetas = np.array(thetas)

    zeta_interp = np.interp(theta_theory, thetas, zetas)

    Ft = -8/15 * np.log(zeta_interp) + 0.9588 # eq 5
    Fr = 2/15 * np.log(zeta_interp) + 0.2526 # eq 6
    Tt = 1/10 * np.log(zeta_interp) + 0.1895 # eq 8
    Tr = -2/5 * np.log(zeta_interp) + 0.3817 # eq 9
    W = data['f_gravity_mag'] # pN = p kg m s^-2
    W /= 1e6 # kg um s^-2

    # below from eqs 4 and 7 with Ff and Fn = 0
    theory_zhao = W * np.sin(theta_theory * np.pi/180) / (6 * np.pi * eta * particle_diameter/2) * (Ft - Fr * Tt / Tr)**(-1)
    ax.plot(theta_theory, theory_zhao, label='theory [Zhao 2002]')

    eleanor = 0.616 * np.sin(theta_theory*np.pi/180) # 0.616 interpreted from graph from Eleanor in slack
    ax.plot(theta_theory, eleanor, label='experiment')

    ax.legend()
    ax.set_xlabel(r'$\theta$ ($^\circ$)')
    ax.set_ylabel('drift x ($\mathrm{\mu m/s}$)')

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1)
    files = common.files_from_argv('particle_linking/data', 'trajs_')
    go(files, ax)
    common.save_fig(fig, f'particle_linking/figures_png/drift_vs_theta_{files[0]}.png')