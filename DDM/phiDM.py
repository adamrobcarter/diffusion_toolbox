import numpy as np
import common
import scipy.optimize
import tqdm
import warnings
from preprocessing.stack_movie import save_array_movie

def calc(stack, pixel_size):
    use_every_nth_frame = 10
    stack = stack[::use_every_nth_frame, :, :]

    stack = stack - stack.mean(axis=0)
    

    print(common.arraysize(stack))
    u_x, u_y, I = common.fourier_2D(stack, axes=(1, 2), spacing=pixel_size)
    # u_xx, u_yy = np.meshgrid(u_x, u_y)
    # for eleanor0.01, I comes out at 26GB! it has dtype complex128
    phi = np.angle(I)
    print(phi.max(), phi.min())
    del I
    d_phi = phi[1:, :, :] - phi[:-1, :, :]
    # d_phi = phi[:, :, :] - phi[0, :, :]
    print(d_phi.max(), d_phi.min())
    # save_array_movie(phi, f'DDM/figures_png/phi.gif')
    del phi
    print(common.arraysize(d_phi))
    d_phi = np.mod(d_phi, 2*np.pi) # because of the subtraction d_phi is in [-2pi, 2pi], so first correct that
    d_phi = np.unwrap(d_phi[:, ::1, ::1], axis=1)
    d_phi = d_phi - d_phi.mean() # we need zero to be in the middle, so the "plane" is centered at the origin
    save_array_movie(d_phi[:10, :, :], f'DDM/figures_png/d_phi.gif')
    print(d_phi[0, :, :].mean(axis=1))

    print(d_phi.sum(axis=(1, 2)))

    d_Rs = np.full((d_phi.shape[0], 2), np.nan)
    
    for t in tqdm.trange(d_phi.shape[0]):

        func = lambda d_R: np.abs( d_phi[t, :, :] - (u_x * d_R[0] + u_y * d_R[1]) ).sum()
        res = scipy.optimize.minimize(func, x0=(0.1, 0.1))

        if res.message != 'Optimization terminated successfully.':
            warnings.warn(res.message)

        if res.x[0] == res.x[1]:
            warnings.warn('x and y were equal')

        print(res.x)

        d_Rs[t, :] = res.x

    return d_Rs

