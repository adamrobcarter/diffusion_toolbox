import numpy as np
import scipy.stats
import multiprocessing
import functools
import common
import tqdm

def intermediate_scattering(log, F_type, crop, num_k_bins, num_iters, d_frames, data, num_frames, max_K, width, height):
    assert not np.isnan(max_K)

    # data is x,y,t
    data[:, 2] -= data[:, 2].min() # convert time to being 0-based

    Fs    = np.full((len(d_frames), num_k_bins), np.nan)
    F_unc = np.full((len(d_frames), num_k_bins), np.nan)
    print(f'F size {common.arraysize(Fs)}')
    ks   = np.full((len(d_frames), num_k_bins+1), np.nan) # +1 because we get the bin edges

    min_K = 2*np.pi/( min(width, height) * crop )

    # first find the particles at each timestep, otherwise we're transferring
    # the whole of data to each process
    particles_at_frame = []
    for frame in range(num_frames):
        # select only particles at the relevent time step
        particles = data[data[:, 2]==frame, :]
        # select only x and y columns
        particles = particles[:, 0:2]
        particles_at_frame.append(particles)

    parallel = True
    if parallel:
        with multiprocessing.Pool(16) as pool:

            internal_bound = functools.partial(intermediate_scattering_for_dframe, log=log, F_type=F_type, crop=crop,
                                num_k_bins=num_k_bins, num_iters=num_iters, d_frames=d_frames, particles_at_frame=particles_at_frame,
                                num_frames=num_frames, min_K=min_K, max_K=max_K, width=width, height=height)

            computation = pool.imap(internal_bound, range(len(d_frames)), chunksize=1)
            results = list(tqdm.tqdm(computation, total=len(d_frames)))
            
            for dframe_i in range(len(d_frames)):
                F_, F_unc_, k_ = results[dframe_i]
                Fs   [dframe_i, :] = F_
                F_unc[dframe_i, :] = F_unc_
                ks   [dframe_i, :] = k_

    else:
        for dframe_i in tqdm.trange(len(d_frames)):
            F_, F_unc_, k_ = intermediate_scattering_for_dframe(dframe_i, log=log, F_type=F_type, crop=crop,
                                num_k_bins=num_k_bins, num_iters=num_iters, d_frames=d_frames, particles_at_frame=particles_at_frame,
                                num_frames=num_frames, min_K=min_K, max_K=max_K, width=width, height=height)
            
            Fs   [dframe_i, :] = F_
            F_unc[dframe_i, :] = F_unc_
            ks   [dframe_i, :] = k_

    nan_fraction = np.isnan(Fs).sum() / Fs.size
    assert nan_fraction < 0.5, f'F nan fraction was {nan_fraction*100:.0f}%'

    # now remove irrelevent data
    # min_useful_K = 2*np.pi/min(width, height)
    # print('exit    ', ks[0, :6])
    # print('exit    ', ks[0, :6] <= min_K)
    # print(min_useful_K, min_K)
    Fs[ks[:, :-1] <= min_K] = np.nan 
    Fs[ks[:,  1:] >  max_K] = np.nan # remove any greater than max_K (limited data here as get only diagonals)
    #  ks[:, 1:] - this is because we wanna test with the left/right bin edge

    # the first k bin is always all nan, so remove it
    assert np.isnan(Fs[:, 0]).sum() == Fs[:, 0].size
    ks    = ks   [:, 1:]
    Fs    = Fs   [:, 1:]
    F_unc = F_unc[:, 1:]

    # print('min_K', )
    
    ks = (ks[:, 1:] + ks[:, :-1])/2 # return bin midpoints not edges

    return Fs,F_unc,ks

def intermediate_scattering_for_dframe(dframe_i, log, F_type, crop, num_k_bins, num_iters, d_frames, particles_at_frame, num_frames, min_K, max_K, width, height):
    d_frame = int(d_frames[dframe_i])
    offset = (num_frames - d_frame - 1) // num_iters
    assert(num_iters * offset + d_frame < num_frames)

    F = np.full((num_iters, num_k_bins), np.nan)
    k = np.full((num_k_bins+1,),           np.nan) # +1 b/c we get the left and right of the final bin
                    
    if F_type == 'F_s':
        func = self_intermediate_scattering_internal
    else:
        func = intermediate_scattering_internal

    for iter_i in range(num_iters):
        frame = iter_i * offset

        particles_t0, particles_t1 = preprocess_scattering(particles_at_frame[frame], particles_at_frame[frame+d_frame], crop=crop)
        width  = width  * crop
        height = height * crop
    
        k_x, k_y, k_bins = get_k_and_bins_for_intermediate_scattering(min_K, max_K, num_k_bins, log_calc=log, log_bins=log)

        k_unbinned, F_unbinned = func(particles_t0, particles_t1, k_x, k_y)
        
        k_, F_ = postprocess_scattering(k_unbinned, F_unbinned, k_bins)
        
        F[iter_i, :] = F_
        if iter_i == 0:
            k = k_
        else:
            assert np.array_equal(k, k_)

    assert np.isnan(F).sum() < F.size

                # print(f"nan S: {np.isnan(F).sum()/F.size:.2f}, nan k: {np.isnan(k).sum()/k.size:.2f}")
                #assert(np.isnan(S).sum()/S.size < 0.5)
                #print(f'min_K={min_K:.3f}, k bin size={k_[1]-k_[0]:.3f}, num bins={num_k_bins}')

    # need nanmean because binned_statistic will return nan if the bin is empty
    return np.mean(F, axis=0), np.std(F, axis=0), k

def preprocess_scattering(particles_t0, particles_t1, crop):
    # first remove any nan particles (particles that exist at t0 but not t1)
    # t0_nans = np.any(np.isnan(particles_t0), axis=1)
    # t1_nans = np.any(np.isnan(particles_t1), axis=1)
    # # nans = t0_nans | t1_nans
    # # print(f'missing particles: {nans.sum()/nans.size*100}%')
    # particles_t0 = particles_t0[~t0_nans, :]
    # particles_t1 = particles_t1[~t1_nans, :]

    # assert np.isnan(particles_t0).sum() == 0
    # assert np.isnan(particles_t1).sum() == 0

    # then do the crop
    # width_thresh  = ( particles_t0[:, 0].max() - particles_t0[:, 0].min() ) * crop
    # height_thresh = ( particles_t0[:, 1].max() - particles_t0[:, 1].min() ) * crop
    # removed_particles_t0 = (particles_t0[:, 0] > width_thresh) | (particles_t0[:, 1] > height_thresh)
    # removed_particles_t1 = (particles_t1[:, 0] > width_thresh) | (particles_t1[:, 1] > height_thresh)
    # removed_particles = removed_particles_t0 | removed_particles_t1
    # particles_t0 = particles_t0[~removed_particles, :]
    # particles_t1 = particles_t1[~removed_particles, :]

    return particles_t0, particles_t1

def postprocess_scattering(k, F, k_bins):
    assert np.isnan(F).sum() == 0, f'S was {np.isnan(F).sum()/F.size*100:.0f}% NaN'

    # print()
    # print('k_bins  ', k_bins[:5])

    F_binned, k_binned, _ = scipy.stats.binned_statistic(k.flatten(), F.flatten(), 'mean', bins=k_bins)
    # binned statistic returns NaN if the bin is empty
    
    assert np.isnan(F_binned).sum() < F_binned.size

    # print('k_binned', k_binned[:5])

    # best to return the middle of the bin not one side
    # k_binned = (k_binned[1:] + k_binned[:-1])/2

    # print('middled ', k_binned[:5])

    return k_binned, F_binned

def get_k_and_bins_for_intermediate_scattering(min_K, max_K, num_k_bins, log_calc, log_bins):

    if log_calc:
        k_x_pos =  np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins, dtype='float64')
        k_x_neg = -np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins, dtype='float64')
        k_x = np.concatenate((k_x_neg, (0,), k_x_pos))
        k_y = np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins, dtype='float64')
        bin_edges = np.concatenate(((0,), k_x_pos))
        # here we invent a bin 0 < k < min_K. Anything in here should be thrown away later
    else:
        # have checked this starting from min_K not -max_K and it does indeed seem to make no difference
        k_x = np.arange(-max_K, max_K, min_K, dtype='float64')
        k_y = np.arange( 0,     max_K, min_K, dtype='float64')

    bins = bin_edges if log_bins else num_k_bins

    assert np.isnan(k_x).sum() == 0
    assert np.isnan(k_y).sum() == 0

    return k_x, k_y, bins

def intermediate_scattering_internal(particles_t0, particles_t1, k_x, k_y):
    # Thorneywork et al 2018 eq (27))

    num_particles = particles_t0.shape[0]
    #print(f"kept {num_particles} of {num_particles_before}")
    
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
    # k_dot_r_mu = k_x * x_mu + k_y * y_mu
    # k_dot_r_nu = k_x * x_nu + k_y * y_nu

    cos_term1 = np.cos(k_dot_r_mu).sum(axis=(0, 1))
    cos_term2 = np.cos(k_dot_r_nu).sum(axis=(0, 1))
    cos_accum = cos_term1 * cos_term2
    del cos_term1, cos_term2
    
    sin_term1 = np.sin(k_dot_r_mu).sum(axis=(0, 1))
    sin_term2 = np.sin(k_dot_r_nu).sum(axis=(0, 1))
    sin_accum = sin_term1 * sin_term2
    del sin_term1, sin_term2
    del k_dot_r_mu, k_dot_r_nu
    
    contrib = 1/num_particles * ( cos_accum + sin_accum )
    k = np.sqrt(k_x**2 + k_y**2)
    # del cos_accum, sin_accum # probably unneeded

    return k, contrib

def self_intermediate_scattering_internal(particles_t0, particles_t1, k_x, k_y):
    # Thorneywork et al 2018 eq (27)

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
    # TODO: do we not need to consider negative n?!!
    # actually I think not because we already consider u -> v and v -> u

    # k_dot_dr = k_x * (x_mu - x_nu)  +  k_y * (y_mu - y_nu)
    k_dot_dr = np.multiply(k_x, x_mu - x_nu, dtype='float64') + np.multiply(k_y, y_mu - y_nu, dtype='float64')

    S = 1/num_particles * np.cos(k_dot_dr).sum(axis=(0))
    del k_dot_dr
    #neg = np.sum(contrib < 0)
    #print(f'{neg/contrib.size:.2f} negative')

    k = np.sqrt(k_x**2 + k_y**2)

    return k, S