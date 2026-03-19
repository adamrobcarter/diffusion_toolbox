import numpy as np
import common
import tqdm

if __name__ == '__main__':
    # for L in np.logspace(np.log10(10), np.log10(1000), num=20):
    # for L in [1000]:
    if True:

        phi = 0.1
        sigma = 0.5
        dt = 256
        D = common.stokes_einstein_D(sigma) / common.S_k_zero(phi)
        print(f'D={D:.3g}um^2/s')
        max_t = 24 * 60 * 60  # 24 hours
        dt = 256
        max_t = 1 * 60 # 1 min
        dt = 0.1
        num_timesteps = int(max_t / dt)
        pixel_size = 0.2

        mult = 4
        L_original = 207
        L_full     = L_original * mult
        num_particles_original = int(L_original**2 * 4 / np.pi * phi / sigma**2)
        num_particles_full     = int(L_full    **2 * 4 / np.pi * phi / sigma**2)
        print(f'num_particles_full = {num_particles_full}')

        rng = np.random.default_rng()

        print('getting data')
        stepsize = np.sqrt( 2 * D * dt )
        print(f'stepsize {stepsize:.3g}')
        print(f'steps array size {num_particles_full} x {num_timesteps} x f32 = {num_particles_full*num_timesteps*4/1e9:.1f}GB')
        steps_x = stepsize * rng.standard_normal(size=(num_particles_full, num_timesteps), dtype=np.float32)
        print('allocated steps_x')
        startpoints_x = rng.uniform(0, L_full,        size=(num_particles_full),             ).astype(np.float32)
        steps_y = stepsize * rng.standard_normal(size=(num_particles_full, num_timesteps), dtype=np.float32)
        print('allocated steps_y')
        startpoints_y = rng.uniform(0, L_full,        size=(num_particles_full),             ).astype(np.float32)

        print('summing')
        x = startpoints_x[:, np.newaxis] + np.cumsum(steps_x, axis=1)
        y = startpoints_y[:, np.newaxis] + np.cumsum(steps_y, axis=1)

        # build filename
        t_string = f"{common.format_time(max_t)}_dt{common.format_time(dt)}"
        filename = f'sim_nointer_{phi}_L{L_original}_t{t_string}_sigma{sigma}_precropL{L_full}'

        # we need to do the modding into box before building the array
        print('modding into box')
        x = x % L_full
        y = y % L_full

        # form the trackpy style array
        trajs = np.full((num_particles_original*num_timesteps, 4), np.nan, dtype=np.float32)
        row = 0
        for i in tqdm.trange(num_particles_full, desc='forming array'):
            for t in range(num_timesteps):
                if x[i, t] >= L_original or y[i, t] >= L_original:
                    continue

                if row >= trajs.shape[0]:
                    break # this is possible if the particles are not quite evenly distributed across the full and original areas

                trajs[row, :] = [x[i, t], y[i, t], t, i]

                row += 1
        # it's also possible (if periodic bcs) that we haven't filled up particles yet
        assert np.isfinite(trajs).sum() / trajs.size > 0.99, f'nan {np.isnan(trajs).sum() / trajs.size:.1%}'
        trajs = trajs[~np.isnan(trajs).any(axis=1)]

        pack_frac = common.calc_pack_frac(trajs, sigma, L_original, L_original, dimension=2)

        # to save the trajectories without periodic BCs, we gotta split the trajectories when they exit
        # we can assume if a particle leaves, it won't come back the otherside the next frame.
        # so let's just do the cropping, and then any non-contiguous trajectories we split

        # some particles might have never been present in the orginal box, so we re-index
        id_map = {}
        next_free_id = 0
        for row in tqdm.trange(trajs.shape[0], desc='reindexing ids'):
            original_id = int(trajs[row, 3])
            if original_id not in id_map:
                id_map[original_id] = next_free_id
                trajs[row, 3] = next_free_id
                next_free_id += 1
            trajs[row, 3] = id_map[original_id]

        # now split trajectories when non-contiguous
        # particles are currently sorted by ID
        id_map = {}
        for row in tqdm.trange(trajs.shape[0], desc='splitting non-contiguous'):
            current_id = trajs[row, 3]
            if current_id in id_map:
                current_id = id_map[current_id]
                trajs[row, 3] = current_id

            if current_id == trajs[row-1, 3]: # same particle as previous row
                if trajs[row, 2] != trajs[row-1, 2] + 1: # not contiguous in time, so make new id
                    new_id = next_free_id
                    id_map[current_id] = new_id
                    trajs[row, 3] = new_id
                    next_free_id += 1

        common.save_data(f'particle_linking/data/trajs_{filename}.npz',
            particles=trajs,
            dimension=2,
            window_size_x=L_original, window_size_y=L_original,
            time_step=dt, particle_diameter=sigma, pack_frac_given=phi,
            simulated_D = D, pack_frac = pack_frac, pixel_size=pixel_size, # fake but used in fkt calculation for determining min_k
        )
        
        common.save_data(f'particle_detection/data/particles_{filename}.npz',
            particles=trajs[:, [0, 1, 2]],
            dimension=2,
            window_size_x=L_original, window_size_y=L_original,
            time_step=dt, particle_diameter=sigma, pack_frac_given=phi,
            simulated_D = D, pack_frac = pack_frac, pixel_size=pixel_size, # fake but used in fkt calculation for determining min_k
        )

        common.save_data(f'particle_detection/data/particles_{filename}.npz',
            particles=trajs[:, (0, 1, 2)],
            dimension=2,
            window_size_x=L_original, window_size_y=L_original,
            time_step=dt, particle_diameter=sigma, pack_frac_given=phi,
            simulated_D = D, pack_frac = pack_frac,
        )