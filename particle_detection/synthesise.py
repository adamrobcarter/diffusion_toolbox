import common
import numpy as np
import tqdm

data = common.load('particle_detection/data/particles_channel.npz')

pixel_size = data['pixel_size']
Lx = data['window_size_x']
Ly = data['window_size_y']
num_pixels_x = int(Lx / pixel_size)
num_pixels_y = int(Ly / pixel_size)
# print('pixels', num_pixels_x, num_pixels_y)
sigma = 2.8
radius = sigma/2
time_step = 1/30
particles = data['particles']

num_timesteps = int(particles[:, 2].max()) + 1

# print('x min max', particles[:, 0].min(), particles[:, 0].max())
# print('y min max', particles[:, 1].min(), particles[:, 1].max())

pixel_x = np.arange(0, num_pixels_x) * pixel_size
pixel_y = np.arange(0, num_pixels_y) * pixel_size
pixel_xx, pixel_yy = np.meshgrid(pixel_x, pixel_y, indexing='ij')

# print('pixel x min max', pixel_xx.min(), pixel_xx.max())
# print('pixel y min max', pixel_yy.min(), pixel_yy.max())

# def single_gaussian_particle(x, y):
#     # i = np.zeros([num_pixels, num_pixels])
#     r2 = (x - pixel_xx)**2 + (y - pixel_yy)**2
#     # i[r2 < 3*sigma] = np.exp(- r2[r2 < 3*sigma] / sigma)
#     # return i
#     i = np.exp(- r2 / radius)
#     assert np.isfinite(i).all()
#     return i

# def single_spherical_particle(x, y):
#     map_to_calc_over_x = np.abs(x - pixel_xx) < radius * 2
#     map_to_calc_over_y = np.abs(y - pixel_yy) < radius * 2
#     map = map_to_calc_over_x & map_to_calc_over_y
#     # print(map.sum(), map.size)

#     r2 = (x - pixel_xx[map])**2 + (y - pixel_yy[map])**2
#     r2[r2 > radius**2] = radius**2 # hack to prevent sqrt of negative number?

#     i = np.zeros_like(pixel_xx)
#     i[map][r2 < radius**2] = np.sqrt(radius**2 - r2[r2 < radius**2])

#     # i = np.sqrt(radius**2 - r2)
#     # i /= radius # normalise so max is 1
#     return i


stack = np.zeros([num_timesteps, num_pixels_x, num_pixels_y])

for row_i in tqdm.trange(particles.shape[0]):
    t = int(particles[row_i, 2])

    x, y = particles[row_i, 0], particles[row_i, 1]
    
    map_x = np.abs(x - pixel_x) < radius * 1.5
    map_y = np.abs(y - pixel_y) < radius * 1.5
    submatrix = np.ix_(map_x, map_y)

    r2 = (x - pixel_xx[submatrix])**2 + (y - pixel_yy[submatrix])**2
    r2[r2 > radius**2] = radius**2 # hack to prevent sqrt of negative number?

    intensity_submatrix = np.sqrt(radius**2 - r2)
    assert intensity_submatrix.sum() > 0
    stack[t][submatrix] += intensity_submatrix

assert stack.max() > 0

common.save_data(f'preprocessing/data/stack_channel_synth.npz',
                stack=stack, pixel_size=pixel_size, time_step=time_step,
                particle_diameter=sigma)
common.save_data(f'preprocessing/data/stack_channel_synth_small.npz',
                stack=stack[::int(num_timesteps/50)], pixel_size=pixel_size, time_step=time_step,
                particle_diameter=sigma)