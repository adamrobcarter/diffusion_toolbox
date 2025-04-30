import numpy as np
import common
import tqdm

data = np.loadtxt('raw_data/ryker/colloid_pos.csv', delimiter=',')

num_timesteps = int(data.shape[0])
num_particles = int((data.shape[1]-1) / 3)
num_rows = num_timesteps * num_particles

particles = np.full((num_rows, 5), np.nan)

row = 0
dt = 1

for frame in tqdm.trange(num_timesteps):
    t = data[frame, 0]
    assert t >= 0, f't = {t}'
    for particle in range(num_particles):
        # print('t/frame', t/frame)
        particles[row, [0, 1, 2]] = data[frame, 1+particle*3:1+(particle+1)*3]
        particles[row, 3] = t
        particles[row, 4] = particle

        row += 1

assert np.isfinite(particles).all()

common.save_data(
    'particle_linking/data/trajs_sim_porous0.npz',
    particles=particles, max_time_hours=particles[:, 3].max()/60/60,
    particle_diameter=5.12*2, dimension=3)