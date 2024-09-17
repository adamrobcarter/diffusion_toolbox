import pickle
import common

with open('raw_data/eleanorlong010.p', 'rb') as f:
    x = pickle.load(f)
    particles = x.to_numpy()

particles[:, 2] -= particles[:, 2].min() # make zero-based
particles[:, 3] -= particles[:, 3].min() # make zero-based

pixel_size = 0.288
particles[:, [0, 1]] *= pixel_size

num_timesteps = particles[:, 2].max() + 1
window_size_x = particles[:, 0].max()
window_size_y = particles[:, 1].max()

common.save_data(f'particle_detection/data/particles_eleanorlong010.npz', particles=particles,
        time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
        window_size_x=window_size_x, window_size_y=window_size_y,
        pack_frac_given=0.105,
        num_timesteps=num_timesteps)
# np.save(f'particle_detection/data/particles_eleanorlong.npy', data_param)
common.save_data(f'particle_linking/data/trajs_eleanorlong010.npz', particles=particles,
        time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
        window_size_x=window_size_x, window_size_y=window_size_y,
        pack_frac_given=0.105,
        num_timesteps=num_timesteps)