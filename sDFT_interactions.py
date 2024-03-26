import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import scipy.special
import scipy.interpolate
import common

num_integration_points = 2000

def trapezoid(x, y):
    return scipy.integrate.trapezoid(x, y)

data_theta = common.load('visualisation/theta_int.npz')
kL_over2_values = data_theta['kL_over2']
thetaInt  = data_theta['theta_int']

def sDFT_interactions(L, t, phi, D0, sigma):
    # returns <N^2(t)> / <N>

    Ntplot = np.zeros_like(t)

    k_values = kL_over2_values / L * 2
    
    S = common.structure_factor_2d_hard_spheres(k_values, phi, sigma)

    # S = np.ones_like(k_values)
    # f_v = 4 / (k_values**4 * L**2) * thetaInt
    f_v = L**2 * thetaInt

    if scalar_input := np.isscalar(t): # hack that allows you to input a scalar t not an array
        t = np.array([t])
        Ntplot = np.array([np.nan])
    
    # for time_index, t in tqdm.tqdm(enumerate(t), total=t.size):
    for time_index, t in enumerate(t):
        # intu = np.zeros_like(kL_values) * units.micrometer
        # calculate for all K values the values inside the integral

        # for k_index, kL in enumerate(kL_values):
        #     c_over_sigma2_value = c_over_sigma2(2 * kL / L * sigmav, phiv)
        #     S = 1 / (1 - (4 * phiv) * c_over_sigma2_value / np.pi) # Hansen & McDonald (3.6.10)
        #     # intu[k_index] = S * thetaInt[k_index] / np.pi**2 * kL * ( 1 - np.exp(-4 * D0 * t / L**2 * kL**2 / S)) / kL**4     * 4 * L
        #     k = k_values[k_index]
        #     f_v = 16 / (k**4 * L**2) * thetaInt[k_index]
        #     intu[k_index] = k / (2 * np.pi)**2 * f_v * S * ( 1 - np.exp(-4 * D0 * t / L**2 * kL**2 / S))

        # intu[k_index] = S * thetaInt[k_index] / np.pi**2 * kL * ( 1 - np.exp(-4 * D0 * t / L**2 * kL**2 / S)) / kL**4     * 4 * L
        integrand = k_values / (2 * np.pi)**2 * f_v * S * ( 1 - np.exp(-1 * D0 * t * (k_values)**2 / S)) # countoscope eq. 4
        # interpolate them and then compute the integral
        # funcu = scipy.interpolate.interp1d(k_values, integrand, kind='cubic')
        # Ntplot[time_index] = scipy.integrate.quad(funcu, np.min(k_values).magnitude, np.max(k_values).magnitude, limit=num_integration_points)[0] * k_values.units * integrand.units # this is the plateau values (or 2/N0 Plateau value...)
        # 
        Ntplot[time_index] = trapezoid(integrand, k_values)

    # Calculate plateau value as well (corresponds to tv = 0 in above formula)
    # integrand = np.zeros_like(kL_values)
    # for k_index, kL in enumerate(kL_values):
    #     c_over_sigma2_value = c_over_sigma2(2 * kL / L * sigmav, phiv)
    #     S = 1.0 / (1 - (4 * phiv) * c_over_sigma2_value / np.pi)
    #     integrand[k_index] = S * thetaInt[k_index] / np.pi ** 2 * kL / kL ** 4

    # funcu = scipy.interpolate.interp1d(kL_values, intu, kind='cubic', fill_value="extrapolate")
    # Pval = scipy.integrate.quad(funcu, 0, np.max(kL_values), limit=num_integration_points)[0] # this is the plateau values (or 2/N0 Plateau value...)
    # N0 = phiv * L**2 / (np.pi * sigmav**2 / 4) # average number of particles in the box
    # plot the msd curve
    # plt.plot(tvalues / (L**2 / (4 * D0)), -Ntplot*2 + 2*Pval)

    return Ntplot if not scalar_input else Ntplot[0]# + 2*Pval

    # plt.legend()
    

if __name__ == '__main__':
    sigma = 2.8 # Particle diameter in um
    D0 = 0.038/2.2 # Bare diffusion coefficient in um^2/s -- short time?
    phi = 0.66 # volume fraction of particles
    L = 4 # values of box size in um
    tvalues = np.logspace(-1, 4, 40) # times is s to calculate the msd

    mean_square_fluctuations = sDFT_interactions(L, tvalues, phi, D0, sigma)

    plt.plot(tvalues, mean_square_fluctuations)
    plt.xlabel('t')
    plt.loglog()
    plt.savefig('figures_png/plateaus2.png')