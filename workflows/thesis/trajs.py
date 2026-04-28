import numpy as np
import common
import tqdm
import matplotlib.pyplot as plt
import workflows.thesis.common

NUM_TO_PLOT = 3
phi = 0.01
sigma = 3
dt = 0.5
D = 0.04*100
# num_timesteps = int(24 * 60 * 60 / dt / 24) # 24 hours
max_t = 50
num_timesteps = int(max_t / dt)


fig, ax = plt.subplots(1, 1, figsize=workflows.thesis.common.figsize_small)

num_particles = 10000

rng = np.random.default_rng(seed=0)

print('getting data')
stepsize = np.sqrt( 2 * D * dt )
print(f'stepsize {stepsize:.3g}')
steps_x = rng.normal(0, stepsize, size=(num_particles, num_timesteps))
startpoints_x = np.zeros(num_particles)
steps_y = rng.normal(0, stepsize, size=(num_particles, num_timesteps))
startpoints_y = np.zeros(num_particles)

print('summing')
x = startpoints_x[:, np.newaxis] + np.cumsum(steps_x, axis=1) - steps_x[:, 0][:, np.newaxis] # subtract the first step so that the MSD starts at 0
y = startpoints_y[:, np.newaxis] + np.cumsum(steps_y, axis=1) - steps_y[:, 0][:, np.newaxis]

for i in range(NUM_TO_PLOT):
    ax.plot(dt*np.arange(num_timesteps), x[i, :])

ax.set_ylabel(r'$x(t)$ ($\mathrm{\mu m}$)')
ax.set_xlabel(r'$t$ ($\mathrm{s}$)')


common.save_fig(fig, f'workflows/thesis/figures/trajs.pdf', hide_metadata=True)
common.save_fig(fig, f'workflows/thesis/figures/trajs.png')