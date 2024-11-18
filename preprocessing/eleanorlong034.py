import numpy as np
import common

print('loading')
all_data = np.loadtxt('raw_data/0.34_EKRM_trajs.dat', delimiter=',', skiprows=1) # Eleanor

print('reshaping')
particles = all_data.reshape((-1, 4))
PIXEL = 0.288
particles[:, [0,1]] *= PIXEL

print('x min max', particles[:, 0].min(), particles[:, 0].max())
print('y min max', particles[:, 1].min(), particles[:, 1].max())

particles[:, 2] -= particles[:, 2].min() # make t zero-based
particles[:, 3] -= particles[:, 3].min() # make ID zero-based
num_timesteps = particles[:, 2].max()
print('num_timesteps', num_timesteps)

window_size_x = particles[:, 0].max()
window_size_y = particles[:, 1].max()

particle_diameter = 3.09

print('edge crop')
EDGE_CROP = 3
size_before = particles.size
particles = common.crop_particles(particles, window_size_x-EDGE_CROP, window_size_y-EDGE_CROP, EDGE_CROP, EDGE_CROP)
window_size_x -= 2*EDGE_CROP
window_size_y -= 2*EDGE_CROP
print(f'edge crop kept {particles.size/size_before:.3f}')

print('saving')
common.save_data(f'particle_detection/data/particles_eleanorlong034.npz', particles=particles,
            time_step=0.5, particle_diameter=particle_diameter, pixel_size=PIXEL,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.342,
            # num_timesteps=num_timesteps
            )
# np.save(f'particle_detection/data/particles_eleanorlong.npy', data_param)
common.save_data(f'particle_linking/data/trajs_eleanorlong034.npz', particles=particles,
         time_step=0.5, particle_diameter=particle_diameter, pixel_size=PIXEL,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.342,
         # num_timesteps=num_timesteps
         )
print('done')