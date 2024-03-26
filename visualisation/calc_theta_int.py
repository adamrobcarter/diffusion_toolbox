import numpy as np
import scipy.integrate
import scipy.special
import scipy.interpolate
import tqdm
import matplotlib.pyplot as plt
import matplotlib.cm
import time

num_integration_points = 2000

def calc(kL_over_2_values):
    theta_int = np.zeros_like(kL_over_2_values)
    for iu, kL_over2 in tqdm.tqdm(enumerate(kL_over_2_values), total=kL_over_2_values.shape[0]):
        # thetafunc= lambda theta: (2*np.sin(kL*np.cos(theta)  )/    np.cos(theta) )**2 * (2*np.sin(kL*np.sin(theta)  )/    np.sin(theta) )**2
        thetafunc  = lambda theta: (2*np.sin(kL_over2*np.cos(theta)  )/(kL_over2*2*np.cos(theta)))**2 * (2*np.sin(kL_over2*np.sin(theta)  )/(kL_over2*2*np.sin(theta)))**2
        # thetafunc2 = lambda theta: (2*np.sin(kL*np.cos(theta)/2)/(kL*np.cos(theta)))**2 * (2*np.sin(kL*np.sin(theta)/2)/(kL*np.sin(theta)))**2
        # thetafunc3 = lambda theta: (2*np.sin(kL*np.cos(theta)  )/(   np.cos(theta)))**2 * (2*np.sin(kL*np.sin(theta)  )/(   np.sin(theta)))**2
        theta_int[iu] = scipy.integrate.quad(thetafunc , 0, 2*np.pi, limit=num_integration_points, epsabs=1e-20)[0]
        # divide by 4

        # theta_int[iu] = scipy.integrate.quad(thetafunc2, 0, 2*np.pi, limit=num_integration_points, epsabs=1e-10)[0] # this is pi/2 periodic so cut in 4...
        # theta_int[iu] = 4 / kL**4 * scipy.integrate.quad(thetafunc3, 0, 2*np.pi, limit=num_integration_points, epsabs=1e-10)[0] # this is pi/2 periodic so cut in 4...
        if iu % 100 == 0:
            theta_th = np.linspace(0, 2*np.pi, 400)
            f = (np.log10(kL_over2) - np.log10(kL_over2values.min())) / (np.log10(kL_over2values.max()) - np.log10(kL_over2values.min()))
            plt.plot(theta_th, thetafunc(theta_th), label=f'$kL/2 = {kL_over2:.4f}$', color=matplotlib.cm.plasma(f))
    return theta_int

#all values of K... hopefully this is broad enough.
# if values at initial times seem off, that might have to do with the need
# to check out a bit more values here...
t0 = time.time()

kL_over2values = np.logspace(-3, 8, 2500)
theta_int = calc(kL_over2values)

# x_theory = np.logspace(2, 6)
# line = lambda x : 1e5 * x**-4

plt.legend(fontsize=7)
plt.semilogy()
plt.savefig('visualisation/figures_png/theta_int_internal.png')
plt.clf()

plt.scatter(kL_over2values, theta_int, s=1)
# plt.plot(x_theory, line(x_theory))
plt.loglog()
plt.savefig('visualisation/figures_png/theta_int.png')

# remove points
# to_remove = theta_int > line(kL_over2values)
# theta_int      = theta_int     [~to_remove]
# kL_over2values = kL_over2values[~to_remove]

np.savez('visualisation/theta_int.npz', theta_int=theta_int, kL_over2=kL_over2values,
         num_integration_points=num_integration_points, computation_time=time.time()-t0)

# kLvalues = np.logspace(-6, 4, 3000)
# np.savez('data/theta_int_big.npz', theta_int=calc(kLvalues), kL=kLvalues, num_integration_points=num_integration_points)