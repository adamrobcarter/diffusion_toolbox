import numpy as np
import scipy.integrate
import scipy.special
import scipy.interpolate
import tqdm

num_integration_points = 2000

def calc(kLvalues):
    theta_int = np.zeros_like(kLvalues)
    for iu, kL_over2 in tqdm.tqdm(enumerate(kLvalues), total=kLvalues.shape[0]):
        # thetafunc= lambda theta: (2*np.sin(kL*np.cos(theta)  )/    np.cos(theta) )**2 * (2*np.sin(kL*np.sin(theta)  )/    np.sin(theta) )**2
        thetafunc  = lambda theta: (2*np.sin(kL_over2*np.cos(theta)  )/(kL_over2*2*np.cos(theta)))**2 * (2*np.sin(kL_over2*np.sin(theta)  )/(kL_over2*2*np.sin(theta)))**2
        # thetafunc2 = lambda theta: (2*np.sin(kL*np.cos(theta)/2)/(kL*np.cos(theta)))**2 * (2*np.sin(kL*np.sin(theta)/2)/(kL*np.sin(theta)))**2
        # thetafunc3 = lambda theta: (2*np.sin(kL*np.cos(theta)  )/(   np.cos(theta)))**2 * (2*np.sin(kL*np.sin(theta)  )/(   np.sin(theta)))**2
        theta_int[iu] = scipy.integrate.quad(thetafunc , 0, 2*np.pi, limit=num_integration_points, epsabs=1e-10)[0]
        # divide by 4

        # theta_int[iu] = scipy.integrate.quad(thetafunc2, 0, 2*np.pi, limit=num_integration_points, epsabs=1e-10)[0] # this is pi/2 periodic so cut in 4...
        # theta_int[iu] = 4 / kL**4 * scipy.integrate.quad(thetafunc3, 0, 2*np.pi, limit=num_integration_points, epsabs=1e-10)[0] # this is pi/2 periodic so cut in 4...
    return theta_int
    # this integral is equal to k^4 L^2 / 16 times f_v(k) as given in the countoscope paper, given K = kL
    # now f_v(k)/L^2

#all values of K... hopefully this is broad enough.
# if values at initial times seem off, that might have to do with the need
# to check out a bit more values here...

kL_over2values = np.logspace(-2.5, 5, 500)
np.savez('box_counting/theta_int.npz', theta_int=calc(kL_over2values), kL_over2=kL_over2values, num_integration_points=num_integration_points)

# kLvalues = np.logspace(-6, 4, 3000)
# np.savez('data/theta_int_big.npz', theta_int=calc(kLvalues), kL=kLvalues, num_integration_points=num_integration_points)