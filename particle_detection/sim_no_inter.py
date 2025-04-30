import numpy as np
import common
import tqdm

L = 640
phi = 0.3
sigma = 3
dt = 0.25
D = 0.04
# num_timesteps = int(24 * 60 * 60 / dt / 24) # 24 hours
num_timesteps = 10000

num_particles = int(L**2 * 4 / np.pi * phi / sigma**2)

rng = np.random.default_rng()

print('getting data')
stepsize = np.sqrt( 2 * D * dt )
print(f'stepsize {stepsize:.3g}')
steps_x = rng.normal(0, stepsize, size=(num_particles, num_timesteps))
startpoints_x = rng.uniform(0, L, size=(num_particles),              )
steps_y = rng.normal(0, stepsize, size=(num_particles, num_timesteps))
startpoints_y = rng.uniform(0, L, size=(num_particles),              )

print('summing')
x = startpoints_x[:, np.newaxis] + np.cumsum(steps_x, axis=1)
y = startpoints_y[:, np.newaxis] + np.cumsum(steps_y, axis=1)

print('modding into box')
x = x % L
y = y % L

trajs = np.full((num_particles*num_timesteps, 3), np.nan, dtype=np.float32)
row = 0

for t in tqdm.trange(num_timesteps, desc='forming array'):
    for i in range(num_particles):
        trajs[row, :] = [x[i, t], y[i, t], t]
        row += 1

# we don't save the particle ID number!
# therefore we have to link to get trajectories later
# this is because if we saved it now, when a particle went across
# the border and came back the other side, it would be on the same
# trajectory which will mess up the MSDs.

phi_str = f'{phi*100:.0f}'.zfill(3)
filename = f'sim_nointer_{phi_str}_L{L}_dt{dt}'
common.save_data(f'particle_detection/data/particles_{filename}.npz',
                 particles=trajs[:, (0, 1, 2)],
                 window_size_x=L, window_size_y=L,
                 time_step=dt, particle_diameter=sigma, pack_frac_given=phi)
# common.save_data(f'particle_linking/data/trajs_{filename}.npz',
#                  particles=trajs,
#                  window_size_x=L, window_size_y=L,
#                  time_step=dt, particle_diameter=sigma, pack_frac_given=phi)