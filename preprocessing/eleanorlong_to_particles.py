import numpy as np

all_data = np.loadtxt('raw_data/0.34_EKRM_trajs.dat', delimiter=',', skiprows=1) # Eleanor
data_param = all_data.reshape((-1, 4))
PIXEL = 0.288
data_param[:, [0,1]] *= PIXEL

num_timesteps = data_param[:, 2].max() - data_param[:, 2].min()

np.savez(f'particle_detection/data/particles_eleanorlong.npz', particles=data_param,
            time_step=0.5, particle_diameter=2.8, pixel_size=PIXEL,
            num_timesteps=num_timesteps)
np.save(f'particle_detection/data/particles_eleanorlong.npy', data_param)
print('done')