import common
import matplotlib.pyplot as plt
import numpy as np
import scipy

if __name__ == '__main__':
    filename = lambda nblobs, theta : f'sim_shear0.0_T0_theta{theta}_RHS_nblobs{nblobs}'
    all_nblobs = [42, 162, 642, 2562]
    all_nblobs = [2562]
    # all_thetas = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]
    all_thetas = [0, 2, 4, 6, 8, 10, 12, 14]

    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(r'$\theta$ ($^{\circ}$)')
    ax.set_ylabel(r'drift x $\mathrm{\mu m/s}$')

    for nblobs in all_nblobs:
        thetas = []
        drifts = []
        zetas = []

        for theta in all_thetas:
        # for file in common.files_from_argv('particle_linking/data', 'trajs_'):
            data = common.load(f'particle_linking/data/trajs_{filename(nblobs, theta)}.npz')
            particle_diameter = data['particle_diameter']   
            drift_xyz = common.find_drift(data['particles'], data.get('dimension', 2))
            z_avg = data['particles'][:, 2].mean()
            z_zeta = z_avg / (particle_diameter/2) - 1

            # print('drift', drift_xyz)
            thetas.append(data['theta'])
            drifts.append(drift_xyz[0]) # um/s
            zetas.append(z_zeta)

        ax.scatter(thetas, drifts, label=f'nblobs={nblobs}')

        ax.legend()
        common.save_fig(fig, f'particle_linking/figures_png/drift_vs_theta_nblobs_zoom_0.png')

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
        theory_zhao = - W * np.sin(theta_theory * np.pi/180) / (6 * np.pi * eta * particle_diameter/2) * (Ft - Fr * Tt / Tr)**(-1)
        ax.plot(theta_theory, theory_zhao, label='Zhao 2002')

        ax.legend()
        common.save_fig(fig, f'particle_linking/figures_png/drift_vs_theta_nblobs_zoom_1.png')

        
    eleanor = -0.616 * np.sin(theta_theory * np.pi/180) # 0.616 interpreted from graph from Eleanor in slack
    ax.plot(theta_theory, eleanor, label='Eleanor exp')

    ax.legend()
    common.save_fig(fig, f'particle_linking/figures_png/drift_vs_theta_nblobs_zoom_2.png')