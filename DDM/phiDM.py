import numpy as np
import common
import scipy.optimize

def calc(stack, pixel_size):
    stack = stack - stack.mean(axis=0)

    u_x, u_y, I = common.fourier_2D(stack, axes=(1, 2), spacing=pixel_size)
    phi = np.angle(I)
    d_phi = phi[1:, :, :] - phi[:-1, :, :]

    print(d_phi.sum(axis=(1, 2)))
    
    for t in range(d_phi.shape[0]):

        func = lambda R: np.abs( np.sum(d_phi[t]) - np.sum(u_x * R[0] + u_y * R[1]) )
        res = scipy.optimize.minimize(func, x0=(0, 0))
        print(res)

        break

