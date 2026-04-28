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

num_particles = 100000

rng = np.random.default_rng(seed=2)

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



for i in range(NUM_TO_PLOT*2):
    ax.plot(dt*np.arange(num_timesteps), x[i, :]**2, color='grey', alpha=0.5)

msd = np.mean(x**2 + y**2, axis=0)
print(msd[:5])
ax.plot(dt*np.arange(num_timesteps)[::3], msd[::3], label='simulation', marker='o', linestyle='none')
ax.plot(dt*np.arange(num_timesteps), 4*D*dt*np.arange(num_timesteps), label='theory', linewidth=2)

ax.set_ylabel(r'$\langle x(t)^2 \rangle$ ($\mathrm{\mu m}^2$)')
ax.set_xlabel(r'$t$ ($\mathrm{s}$)')
ax.legend()

common.save_fig(fig, f'workflows/thesis/figures/msd.pdf', hide_metadata=True)
common.save_fig(fig, f'workflows/thesis/figures/msd.png')