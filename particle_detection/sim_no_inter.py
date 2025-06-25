import numpy as np
import common
import tqdm

for L in np.logspace(np.log10(10), np.log10(1000), num=20):
    L = int(L)

    phi = 0.1
    sigma = 3
    dt = 500
    D = 0.04
    # num_timesteps = int(24 * 60 * 60 / dt / 24) # 24 hours
    max_t = 1e7
    t_string = '1e7'
    num_timesteps = int(max_t / dt)

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

    # save the unwrapped trajectories now before wrapping
    phi_str = f'{phi*100:.0f}'.zfill(3)
    filename = f'sim_nointer_{phi_str}_L{L}_t{t_string}'

    trajs = np.full((num_particles*num_timesteps, 4), np.nan, dtype=np.float32)
    row = 0
    for t in tqdm.trange(num_timesteps, desc='forming array'):
        for i in range(num_particles):
            trajs[row, :] = [x[i, t], y[i, t], t, i]
            row += 1
    common.save_data(f'particle_linking/data/trajs_{filename}_unwrap.npz',
        particles=trajs,
        dimension=2,
        window_size_x=L, window_size_y=L,
        time_step=dt, particle_diameter=sigma, pack_frac_given=phi
    )
    common.save_data(f'particle_detection/data/particles_{filename}_unwrap.npz',
        particles=trajs[:, [0, 1, 2]],
        dimension=2,
        window_size_x=L, window_size_y=L,
        time_step=dt, particle_diameter=sigma, pack_frac_given=phi
    )

    print('modding into box')
    x = x % L
    y = y % L

    particles = np.full((num_particles*num_timesteps, 3), np.nan, dtype=np.float32)
    row = 0
    for t in tqdm.trange(num_timesteps, desc='forming array'):
        for i in range(num_particles):
            particles[row, :] = [x[i, t], y[i, t], t]
            row += 1

    # we don't save the particle ID number!
    # therefore we have to link to get trajectories later
    # this is because if we saved it now, when a particle went across
    # the border and came back the other side, it would be on the same
    # trajectory which will mess up the MSDs.

    common.save_data(f'particle_detection/data/particles_{filename}.npz',
                    particles=particles[:, (0, 1, 2)],
                    dimension=2,
                    window_size_x=L, window_size_y=L,
                    time_step=dt, particle_diameter=sigma, pack_frac_given=phi)