import numpy as np
import common

print('loading')
all_data = np.loadtxt('raw_data/0.34_EKRM_trajs.dat', delimiter=',', skiprows=1) # Eleanor

print('reshaping')
data_param = all_data.reshape((-1, 4))
PIXEL = 0.288
data_param[:, [0,1]] *= PIXEL

print('x min max', data_param[:, 0].min(), data_param[:, 0].max())
print('y min max', data_param[:, 1].min(), data_param[:, 1].max())

data_param[:, 2] -= data_param[:, 2].min() # make t zero-based
data_param[:, 3] -= data_param[:, 3].min() # make ID zero-based
num_timesteps = data_param[:, 2].max()
print('num_timesteps', num_timesteps)

window_size_x = data_param[:, 0].max()
window_size_y = data_param[:, 1].max()

print('saving')
common.save_data(f'particle_detection/data/particles_eleanorlong034.npz', particles=data_param,
            time_step=0.5, particle_diameter=2.8, pixel_size=PIXEL,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.342,
            num_timesteps=num_timesteps)
# np.save(f'particle_detection/data/particles_eleanorlong.npy', data_param)
common.save_data(f'particle_linking/data/trajs_eleanorlong034.npz', particles=data_param,
         time_step=0.5, particle_diameter=2.8, pixel_size=PIXEL,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.342,
         num_timesteps=num_timesteps)
print('done')