import pickle
import common

with open('raw_data/eleanor/eleanorlong001.p', 'rb') as f:
    x = pickle.load(f)
    particles = x.to_numpy()

particles[:, 2] -= particles[:, 2].min() # make zero-based
particles[:, 3] -= particles[:, 3].min() # make zero-based

pixel_size = 0.288
particles[:, [0, 1]] *= pixel_size

window_size_x = particles[:, 0].max()
window_size_y = particles[:, 1].max()

particle_diameter = 3.09

EDGE_CROP = 3
particles = common.crop_particles(particles, window_size_x-EDGE_CROP, window_size_y-EDGE_CROP, EDGE_CROP, EDGE_CROP)
window_size_x -= 2*EDGE_CROP
window_size_y -= 2*EDGE_CROP

common.save_data(f'particle_detection/data/particles_eleanorlong001.npz', particles=particles,
        time_step=0.5, particle_diameter=particle_diameter, pixel_size=pixel_size,
        window_size_x=window_size_x, window_size_y=window_size_y,
        pack_frac_given=0.015,
)
# np.save(f'particle_detection/data/particles_eleanorlong.npy', data_param)
common.save_data(f'particle_linking/data/trajs_eleanorlong001.npz', particles=particles,
        time_step=0.5, particle_diameter=particle_diameter, pixel_size=pixel_size,
        window_size_x=window_size_x, window_size_y=window_size_y,
        pack_frac_given=0.015,
)