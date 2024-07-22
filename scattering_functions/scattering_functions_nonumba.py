import numpy as np
import scipy.stats
import multiprocessing
import functools
import common
import tqdm
import warnings, time

def intermediate_scattering(log, F_type, num_k_bins, max_time_origins, d_frames, particles, max_K, min_K):
    assert not np.isnan(max_K)

    # data is x,y,t
    particles[:, 2] -= particles[:, 2].min() # convert time to being 0-based

    Fs    = np.full((len(d_frames), num_k_bins), np.nan)
    F_unc = np.full((len(d_frames), num_k_bins), np.nan)
    # print(f'F size {common.arraysize(Fs)}')
    ks   = np.full((len(d_frames), num_k_bins+1), np.nan) # +1 because we get the bin edges

    # first find the particles at each timestep, otherwise we're transferring
    # the whole of data to each process
    print('finding particles at each timestep')

    num_timesteps = int(particles[:, 2].max()) + 1
    # num_timesteps = int(particles[:, 2].max())

    # first find max number of particles at any one timestep
    num_particles_at_frame = np.bincount(particles[:, 2].astype('int'))
    max_particles_at_frame = num_particles_at_frame.max()
    print('max particles at any one frame', max_particles_at_frame)

    if F_type == 'F_s':
        # for Fself, we need the IDs, so we provide a list where the nth element is the nth particle
        assert particles.shape[1] == 4, 'for self intermediate scattering, you should provide rows of x,y,t,#'

        # make sure IDs are 0-based
        particles[:, 3] -= particles[:, 3].min()

        # there's probably a quicker way to do the below loop - see num_particles_at_frame
        raise Exception('you need to convert this so that particles_at_frame is an ndarray')

        num_particles = int(particles[:, 3].max()) + 1
        assert num_particles > 0
        # particles_at_frame = [np.full((num_particles, 2), np.nan) for _ in range(num_timesteps)]
        particles_at_frame = [np.full((num_particles, 4), np.nan) for _ in range(num_timesteps)]
        for i in tqdm.trange(len(particles)):
            t = int(particles[i, 2])
            id = int(particles[i, 3])
            particles_at_frame[t][id, :] = particles[i, :]

        for t in tqdm.trange(num_timesteps):
            assert np.isnan(particles_at_frame[t]).sum() < particles_at_frame[t].size, f'particles_at_frame[{t}] was all nan'

    else:
        # for F (not self), we don't need the ID, so we just provide a list of particles
        # some datasets may already have the ID in column 4, so we only select the first 3 columns
        particles = particles[:, [0, 1, 2]]

        # the below is a heavily-optimised method for turning the array of particles
        # into an array that is (num timesteps) x (max particles per timestep) x 2
        # we add extra nan rows such that each timestep has the same number of rows
        num_extra_rows = (max_particles_at_frame - num_particles_at_frame).sum()
        extra_rows = np.full((num_extra_rows, 3), np.nan)
        num_rows_added = 0
        for frame in range(num_timesteps):
            for i in range(max_particles_at_frame-num_particles_at_frame[frame]):
                extra_rows[num_rows_added, 2] = frame
                num_rows_added += 1
        all_rows = np.concatenate((particles, extra_rows), axis=0)
        # then by sorting and reshaping, we can get the structure we want
        all_rows_sorted = all_rows[all_rows[:, 2].argsort()]
        particles_at_frame = all_rows_sorted.reshape((num_timesteps, max_particles_at_frame, 3))
        # now remove the time column, leaving just x and y
        particles_at_frame = particles_at_frame[:, :, [0, 1]]
        # remove the nan rows too
        # @TODO:
        # nans = np.isnan(particles_at_frame[:, :, 0])
        # print(nans.shape)
        # print(particles_at_frame.shape)
        # particles_at_frame = particles_at_frame[~nans, :]
        # print(particles_at_frame.shape)

    del particles

    use_every_nth_frame = max(int(num_timesteps / max_time_origins), 1)

    if use_every_nth_frame > 1:
        warnings.warn(f'Using every {use_every_nth_frame}th frame as a time origin. Eventually you may want to use every frame')

    print('beginning computation')
    print('particles_at_frame:', common.arraysize(particles_at_frame))

    parallel = False
    if parallel:
        # cores = 8
        # print('total_RAM', common.arraysize(particles_at_frame, cores))
        # if cores > 16:
        #     warnings.warn(f'using {cores} cores')
        # with multiprocessing.Pool(cores) as pool:

        #     internal_bound = functools.partial(intermediate_scattering_for_dframe, log=log, F_type=F_type,
        #                         num_k_bins=num_k_bins, use_every_nth_frame=use_every_nth_frame, d_frames=d_frames, particles_at_frame=particles_at_frame,
        #                         num_frames=num_timesteps, min_K=min_K, max_K=max_K)

        #     computation = pool.imap(internal_bound, range(len(d_frames)), chunksize=1)
        #     results = list(tqdm.tqdm(computation, total=len(d_frames)))
            
        #     for dframe_i in range(len(d_frames)):
        #         F_, F_unc_, k_ = results[dframe_i]
        #         Fs   [dframe_i, :] = F_
        #         F_unc[dframe_i, :] = F_unc_
        #         ks   [dframe_i, :] = k_
                
        # nan_fraction = common.nanfrac(Fs)
        # assert nan_fraction < 0.5, f'F nan fraction was {nan_fraction*100:.0f}%'
        pass
    else:
        for dframe_i in tqdm.trange(len(d_frames)):
            F_, F_unc_, k_ = intermediate_scattering_for_dframe(dframe_i, log=log, F_type=F_type,
                                num_k_bins=num_k_bins, use_every_nth_frame=use_every_nth_frame, d_frames=d_frames, particles_at_frame=particles_at_frame,
                                num_frames=num_timesteps, min_K=min_K, max_K=max_K)
            
            Fs   [dframe_i, :] = F_
            F_unc[dframe_i, :] = F_unc_
            ks   [dframe_i, :] = k_

            nan_fraction = common.nanfrac(F_)
            assert nan_fraction < 0.5, f'F nan fraction was {nan_fraction*100:.0f}% for d_frame={d_frames[dframe_i]}'


    # now remove irrelevent data
    # min_useful_K = 2*np.pi/min(width, height)
    Fs[ks[:, :-1] <= min_K] = np.nan 
    Fs[ks[:,  1:] >  max_K] = np.nan # remove any greater than max_K (limited data here as get only diagonals)
    #  ks[:, 1:] - this is because we wanna test with the left/right bin edge

    # the first k bin is always all nan, so remove it
    assert np.isnan(Fs[:, 0]).sum() == Fs[:, 0].size
    ks    = ks   [:, 1:]
    Fs    = Fs   [:, 1:]
    F_unc = F_unc[:, 1:]
    
    ks = (ks[:, 1:] + ks[:, :-1])/2 # return bin midpoints not edges

    return Fs,F_unc,ks

def intermediate_scattering_for_dframe(dframe_i, log, F_type, num_k_bins, use_every_nth_frame, d_frames, particles_at_frame, num_frames, min_K, max_K):
    d_frame = int(d_frames[dframe_i])
    # offset = (num_frames - d_frame - 1) // num_iters
    # assert(num_iters * offset + d_frame < num_frames)

    # use_every_nth_frame = 20
    # if use_every_nth_frame != 1:
    #     pass

    assert num_frames > d_frame, f'd_frame={d_frame}, num_frames={num_frames}'

    frames_to_use = range(0, num_frames-d_frame-1, use_every_nth_frame)

    assert max(frames_to_use) + d_frame < num_frames
    
    num_used_frames = len(frames_to_use)
    # print(f'at d_frame={d_frame} using {num_used_frames} frames')

    F = np.full((num_used_frames, num_k_bins), np.nan)
    k = np.full((num_k_bins+1,),           np.nan) # +1 b/c we get the left and right of the final bin
    # F2 = np.full((num_used_frames, num_k_bins), np.nan)
    # k2 = np.full((num_k_bins+1,),           np.nan) # +1 b/c we get the left and right of the final bin
                    
    if F_type == 'F_s':
        func = self_intermediate_scattering_internal
    else:
        func = intermediate_scattering_internal

    parallel = True

    if not parallel:
        print('note: u are not paralllel')
        for frame_index in range(num_used_frames):
            frame = int(frames_to_use[frame_index])

            particles = (particles_at_frame[frame, :, :], particles_at_frame[frame+d_frame, :, :])
            k_, F_ = intermediate_scattering_internal_internal(min_K, max_K, num_k_bins, log, func, particles)
            
            F[frame_index, :] = F_
            if frame_index == 0:
                k = k_
            else:
                assert np.array_equal(k, k_)

    else:   
        cores = 32
        if cores > 16:
            warnings.warn(f'using {cores} cores')

        with multiprocessing.Pool(cores) as pool:
            bound = functools.partial(intermediate_scattering_internal_internal,
                                      min_K, max_K, num_k_bins, log, func)
            
            particles = []
            for frame_index in range(num_used_frames):
                frame = int(frames_to_use[frame_index])
                particles.append((particles_at_frame[frame, :, :], particles_at_frame[frame+d_frame, :, :]))

            results = pool.map(bound, particles, chunksize=1)

            # results is now (num used frames) x 2 x (len of slice)
            for i, result in enumerate(results):
                F[i, :] = result[1]
                if i == 0:
                    k = result[0]
                else:
                    assert np.array_equal(k, result[0])

                # nan_fraction = common.nanfrac(result[1])
                # # print('nanfrac', nan_fraction)
                # assert nan_fraction < 0.1, f'F nan fraction was {nan_fraction*100:.0f}%'

    # assert np.isnan(F).sum() < F.size

    # print(f"nan S: {np.isnan(F).sum()/F.size:.2f}, nan k: {np.isnan(k).sum()/k.size:.2f}")
                #assert(np.isnan(S).sum()/S.size < 0.5)
                #print(f'min_K={min_K:.3f}, k bin size={k_[1]-k_[0]:.3f}, num bins={num_k_bins}')

    # need nanmean because binned_statistic will return nan if the bin is empty
    return np.nanmean(F, axis=0), np.nanstd(F, axis=0)/np.sqrt(num_used_frames), k

def intermediate_scattering_internal_internal(min_K, max_K, num_k_bins, log, func, particles):
    
    particles_t0, particles_t1 = preprocess_scattering(particles[0], particles[1])

    k_x, k_y, k_bins = get_k_and_bins_for_intermediate_scattering(min_K, max_K, num_k_bins, log_calc=log, log_bins=log)

    k_unbinned, F_unbinned = func(particles_t0, particles_t1, k_x, k_y)
    
    assert np.isnan(F_unbinned).sum() == 0, f'F was {np.isnan(F_unbinned).sum()/F_unbinned.size*100:.0f}% NaN'
    
    k_, F_ = postprocess_scattering(k_unbinned, F_unbinned, k_bins)
    return k_, F_

def preprocess_scattering(particles_t0, particles_t1):
    # first remove any nan particles
    t0_nans = np.any(np.isnan(particles_t0), axis=1)
    t1_nans = np.any(np.isnan(particles_t1), axis=1)
    # # nans = t0_nans | t1_nans
    # # print(f'missing particles: {nans.sum()/nans.size*100}%')
    particles_t0 = particles_t0[~t0_nans, :]
    particles_t1 = particles_t1[~t1_nans, :]

    # assert np.isnan(particles_t0).sum() == 0
    # assert np.isnan(particles_t1).sum() == 0

    return particles_t0, particles_t1

def postprocess_scattering(k, F, k_bins):
    # print()

    F_binned, k_binned, _ = scipy.stats.binned_statistic(k.flatten(), F.flatten(), 'mean', bins=k_bins)
    # binned statistic returns NaN if the bin is empty
    
    assert np.isnan(F_binned).sum() < F_binned.size

    # print('k_binned', k_binned[:5])

    # best to return the middle of the bin not one side
    # k_binned = (k_binned[1:] + k_binned[:-1])/2

    # print('middled ', k_binned[:5])

    return k_binned, F_binned

def get_k_and_bins_for_intermediate_scattering(min_K, max_K, num_k_bins, log_calc, log_bins):

    sym = False

    if log_calc:
        k_x_pos =  np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins, dtype='float64')
        k_x_neg = -np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins, dtype='float64')
        bin_edges = np.concatenate(((0,), k_x_pos))

        k_x = np.concatenate((k_x_neg, (0,), k_x_pos))
        if sym:
            k_y = np.copy(k_x)
        else:
            k_y = np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins, dtype='float64')
        # k_y = np.copy(k_x)
        # here we invent a bin 0 < k < min_K. Anything in here should be thrown away later
    else:
        # have checked this starting from min_K not -max_K and it does indeed seem to make no difference
        k_x_pos = np.linspace(min_K, max_K, num_k_bins, dtype='float64')
        k_x_neg = -np.copy(k_x_pos)
        k_x = np.concatenate((k_x_neg, (0,), k_x_pos))
        # ^^ this is because if we do arange(-max_K, max_K, min_K), if min_K is not a divisor of max_K, k will not
        # be symetrical. I don't know if that's an issue right now, but it could be
        if sym:
            k_y = np.copy(k_x)
        else:
            k_y = np.concatenate(((0,), np.linspace(min_K, max_K, num_k_bins, dtype='float64')))

    bins = bin_edges if log_bins else num_k_bins

    assert np.isnan(k_x).sum() == 0
    assert np.isnan(k_y).sum() == 0
    
    assert np.isclose(k_x.mean(), 0), f'k_x.mean() = {k_x.mean()}'

    return k_x, k_y, bins

def intermediate_scattering_internal(particles_t0, particles_t1, k_x, k_y):
    # Thorneywork et al 2018 eq (27))
    
    particle_t0_x = particles_t0[:, 0]
    particle_t0_y = particles_t0[:, 1]
    particle_t1_x = particles_t1[:, 0]
    particle_t1_y = particles_t1[:, 1]

    # dimensions are
    # mu x nu x kx x ky
    x_mu = particle_t0_x[:, np.newaxis, np.newaxis, np.newaxis]
    y_mu = particle_t0_y[:, np.newaxis, np.newaxis, np.newaxis]
    x_nu = particle_t1_x[np.newaxis, :, np.newaxis, np.newaxis]
    y_nu = particle_t1_y[np.newaxis, :, np.newaxis, np.newaxis]

    k_x = k_x[np.newaxis, np.newaxis, :, np.newaxis]
    k_y = k_y[np.newaxis, np.newaxis, np.newaxis, :]
    
    # TODO: do we not need to consider negative n?!!
    # actually I think not because we already consider u -> v and v -> u

    k_dot_r_mu = np.multiply(k_x, x_mu, dtype='float64') + np.multiply(k_y, y_mu, dtype='float64')
    k_dot_r_nu = np.multiply(k_x, x_nu, dtype='float64') + np.multiply(k_y, y_nu, dtype='float64')
    # print(k_x.mean(), k_y.mean(), x_mu.mean(), x_nu.mean(), y_mu.mean(), y_nu.mean(), k_dot_r_mu.mean(), k_dot_r_nu.mean())
    # k_dot_r_mu = k_x * x_mu + k_y * y_mu
    # k_dot_r_nu = k_x * x_nu + k_y * y_nu

    # print(np.cos(k_dot_r_mu).sum(axis=(0, 1)).shape, np.cos(k_dot_r_mu).sum(axis=(0)).shape)
    cos_term1 = np.cos(k_dot_r_mu).sum(axis=0) # sum over mu
    cos_term2 = np.cos(k_dot_r_nu).sum(axis=1) # sum over nu
    cos_accum = cos_term1 * cos_term2
    del cos_term1, cos_term2
    
    sin_term1 = np.sin(k_dot_r_mu).sum(axis=0)
    sin_term2 = np.sin(k_dot_r_nu).sum(axis=1)
    sin_accum = sin_term1 * sin_term2
    del sin_term1, sin_term2
    del k_dot_r_mu, k_dot_r_nu
    
    num_particles = (particles_t0.shape[0] + particles_t1.shape[0]) / 2
    if num_particles == 0:
        warnings.warn('found no particles in either timestep')
        contrib = np.zeros_like(cos_accum)
    else:
        contrib = 1/num_particles * ( cos_accum + sin_accum )
    k = np.sqrt(k_x**2 + k_y**2)
    # del cos_accum, sin_accum # probably unneeded

    return k, contrib

def distinct_intermediate_scattering_internal(particles_t0, particles_t1, k_x, k_y):
    
    particle_t0_x = particles_t0[:, 0]
    particle_t0_y = particles_t0[:, 1]
    particle_t1_x = particles_t1[:, 0]
    particle_t1_y = particles_t1[:, 1]

    # dimensions are
    # mu x nu x kx x ky
    x_mu = particle_t0_x[:, np.newaxis, np.newaxis, np.newaxis]
    y_mu = particle_t0_y[:, np.newaxis, np.newaxis, np.newaxis]
    x_nu = particle_t1_x[np.newaxis, :, np.newaxis, np.newaxis]
    y_nu = particle_t1_y[np.newaxis, :, np.newaxis, np.newaxis]

    k_x = k_x[np.newaxis, np.newaxis, :, np.newaxis]
    k_y = k_y[np.newaxis, np.newaxis, np.newaxis, :]
    
    # TODO: do we not need to consider negative n?!!
    # actually I think not because we already consider u -> v and v -> u

    k_dot_r_mu = np.multiply(k_x, x_mu, dtype='float64') + np.multiply(k_y, y_mu, dtype='float64')
    k_dot_r_nu = np.multiply(k_x, x_nu, dtype='float64') + np.multiply(k_y, y_nu, dtype='float64')
    # print(k_x.mean(), k_y.mean(), x_mu.mean(), x_nu.mean(), y_mu.mean(), y_nu.mean(), k_dot_r_mu.mean(), k_dot_r_nu.mean())
    # k_dot_r_mu = k_x * x_mu + k_y * y_mu
    # k_dot_r_nu = k_x * x_nu + k_y * y_nu

    # print(np.cos(k_dot_r_mu).sum(axis=(0, 1)).shape, np.cos(k_dot_r_mu).sum(axis=(0)).shape)
    cos_term1 = np.cos(k_dot_r_mu).sum(axis=0) # sum over mu
    cos_term2 = np.cos(k_dot_r_nu).sum(axis=1) # sum over nu
    cos_accum = cos_term1 * cos_term2
    del cos_term1, cos_term2
    
    sin_term1 = np.sin(k_dot_r_mu).sum(axis=0)
    sin_term2 = np.sin(k_dot_r_nu).sum(axis=1)
    sin_accum = sin_term1 * sin_term2
    del sin_term1, sin_term2
    del k_dot_r_mu, k_dot_r_nu
    
    num_particles = (particles_t0.shape[0] + particles_t1.shape[0]) / 2
    if num_particles == 0:
        warnings.warn('found no particles in either timestep')
        contrib = np.zeros_like(cos_accum)
    else:
        contrib = 1/num_particles * ( cos_accum + sin_accum )
    k = np.sqrt(k_x**2 + k_y**2)
    # del cos_accum, sin_accum # probably unneeded

    return k, contrib

def self_intermediate_scattering_internal(particles_t0, particles_t1, k_x, k_y):
    # remove particles that were nan in both sets
    nans = np.isnan(particles_t0[:, 0]) | np.isnan(particles_t1[:, 0])
    
    particles_t0 = particles_t0[~nans, :]
    particles_t1 = particles_t1[~nans, :]
    
    num_particles = particles_t0.shape[0]
    #print(f"kept {num_particles} of {num_particles_before}")
    
    particle_t0_x = particles_t0[:, 0]
    particle_t0_y = particles_t0[:, 1]
    particle_t1_x = particles_t1[:, 0]
    particle_t1_y = particles_t1[:, 1]

    # dimensions are
    # mu x kx x ky
    x_mu = particle_t0_x[:, np.newaxis, np.newaxis]
    y_mu = particle_t0_y[:, np.newaxis, np.newaxis]
    x_nu = particle_t1_x[:, np.newaxis, np.newaxis]
    y_nu = particle_t1_y[:, np.newaxis, np.newaxis]

    k_x = k_x[np.newaxis, :, np.newaxis]
    k_y = k_y[np.newaxis, np.newaxis, :]

    # k_dot_dr = k_x * (x_mu - x_nu)  +  k_y * (y_mu - y_nu)
    k_dot_dr = np.multiply(k_x, x_mu - x_nu, dtype='float64') + np.multiply(k_y, y_mu - y_nu, dtype='float64')

    S = 1/num_particles * np.cos(k_dot_dr).sum(axis=(0))

    del k_dot_dr
    #neg = np.sum(contrib < 0)
    #print(f'{neg/contrib.size:.2f} negative')

    k = np.sqrt(k_x**2 + k_y**2)

    return k, S