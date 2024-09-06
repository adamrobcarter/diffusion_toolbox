import pickle
import common

with open('raw_data/eleanorlong010.p', 'rb') as f:
    x = pickle.load(f)
    particles = x.to_numpy()

particles[:, 2] -= particles[:, 2].min() # make zero-based
particles[:, 3] -= particles[:, 3].min() # make zero-based
num_timesteps = particles[:, 2].max() + 1
pixel_size = 0.288

common.save_data(f'particle_detection/data/particles_eleanorlong010.npz', particles=particles,
        time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
        num_timesteps=num_timesteps)
# np.save(f'particle_detection/data/particles_eleanorlong.npy', data_param)
common.save_data(f'particle_linking/data/trajs_eleanorlong010.npz', particles=particles,
        time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
        num_timesteps=num_timesteps)