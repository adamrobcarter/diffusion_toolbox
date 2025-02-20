import numpy as np
import common
import scipy.stats
import tqdm
import numba
import warnings
import time

@numba.njit
def do_binning(slice_length, F_D_sq_this, u_flat, u_bins):
    F_D_sq = np.zeros((slice_length, len(u_bins)-1))
    for t_index in numba.prange(slice_length):
        # assert u.shape == F_D_sq_this[t_index, :, :].shape
        F_D_sq[t_index, :], _, _ = common.numba_binned_statistic(u_flat, F_D_sq_this[t_index, :, :].flatten(), bins=u_bins)
    return F_D_sq


def calc(stack, pixel_size, time_step, num_k_bins, callback=lambda i, k, F_D_sq, F_D_sq_unc, t, F_D_sq_all, time_origins, F_D_sq_noradial, k_x, k_y : None):
    # num_k_bins is only used if the smallest image dimension > 150
    # otherwise the bins will be linearly (not log) spaced between min_K and max_K
     
    I = stack # we use Cerbino's notation
    assert I.ndim == 3 # catch incase any 2-channel file has got through
    assert common.nanfrac(I) == 0

    use_every_nth_frame = 1
    while stack.size/use_every_nth_frame > 4e8:
        use_every_nth_frame *= 2
    if use_every_nth_frame > 1:
        warnings.warn(f'using every {use_every_nth_frame}th frame')
    # use_every_nth_frame = max(int(stack.size / 1e8), 1)
    # print(f'automatic: use every {use_every_nth_frame}th frame')
    # use_every_nth_frame = 1 # 10 is a good number for Alice, 1 for Marine]
    
    # use_every_nth_frame = max(num_timesteps / max_time_origins, 1) this is what we use in scattering_functions btw

    time_origins = range(0, I.shape[0], use_every_nth_frame)
    print(f'num time origins: {len(time_origins)}')

    used_times = common.exponential_integers(1, I.shape[0]-1) - 1

    smallest_image_dim = min(stack.shape[1], stack.shape[2])

    assert num_k_bins <= smallest_image_dim

    min_K = 2*np.pi/( smallest_image_dim * pixel_size )

    if smallest_image_dim < 150:
        max_K = 2*np.pi/ pixel_size
        num_k_bins = 20
        k_bins = np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins)
    else:
        max_K = 10
        k_bins = np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins)
    # assert k_bins[1] - k_bins[0] > min_K/2, f'you probably have num_k_bins too high or max_k too low (k_spacing[0]={k_bins[1] - k_bins[0]:.3f}, min_K={min_K:.3f})'
    u_bins = k_bins / (2*np.pi)
    num_u_bins = u_bins.size - 1 # minus 1 because we specify right and left of final bin
    
    F_D_sq = np.full((len(time_origins), used_times.size, num_u_bins), np.nan)
    print('F_D_sq array size', common.arraysize(F_D_sq))

    # a = 0
    # b = 0


    # print('doing big fourier')
    # # instead of taking the intensity differences and fourier transforming them, we could fourier transform the intensities
    # # and then take the differences. This is okay cause the fourier transform is linear. However we don't do that cause the 
    # # speed bottleneck is binned_statistic, not the fourier transform
    # u_x, u_y, I_fourier = common.fourier_2D(I, spacing=pixel_size, axes=(1, 2))
    # u = np.sqrt(u_x**2 + u_y**2)
    # print('done')


    F_D_sq_noradial_sum = np.zeros([used_times.size, I.shape[1], I.shape[2]])
    noradial_sum_i = 0
    # again, would be easier for this to be (time origins) x (used times) x Ix x Iy, but then the array would be too big
    # instead we do averaging as we go

    for time_origin_index, time_origin in tqdm.tqdm(enumerate(time_origins), total=len(time_origins)):

        time_indexes = time_origin+used_times
        time_indexes = time_indexes[time_indexes < I.shape[0]] # remove any that are longer than the data
        slice_length = time_indexes.size

        D = I[time_indexes, :, :] - I[time_origin, :, :] # need -1 cause used_times[0]==1

        u_x, u_y, F_D = common.fourier_2D(D, spacing=pixel_size, axes=(1, 2))
        # print('u spacing', u_x[0, 1] - u_x[0, 0], u_y[1, 0] - u_y[0, 0])
        # print('u max', u_x.max(), u_y.max(), 'u_bins max', u_bins.max())
        # print(u_x.shape, u_y.shape, F_D.shape)
        # print('u bins', u_bins)

        # F_D = I_fourier[time_indexes, :, :] - I_fourier[time_origin, :, :]

        warnings.warn('see static_fourier.py, should you be using fftshift?')
    
        F_D_sq_this = np.abs(F_D)**2
        # print('how real:', np.nanmean(np.imag(F_D))/np.nanmean(np.real(F_D)), np.nanmean(np.imag(F_D)), np.nanmean(np.real(F_D)))
    
        u = np.sqrt(u_x**2 + u_y**2)

        F_D_sq_noradial_sum[:F_D_sq_this.shape[0], :, :] += F_D_sq_this
        noradial_sum_i      += 1
        
        # for i in range(0, 10):
        #     print()
        #     if i == 0:
        #         lower = 0
        #     else:
        #         lower = u_bins[i-1]
        #     upper = u_bins[i]
        #     print(lower, upper)
        #     ins = ( lower < u ) & ( u < upper )
        #     print(u_x[ins], u_y[ins])
        # break

        # radial average
        # this is very surely the time bottleneck
        u_flat = u.flatten()
        for t_index in range(slice_length):
            F_D_sq[time_origin_index, t_index, :], _, _ = scipy.stats.binned_statistic(u_flat, F_D_sq_this[t_index, :, :].flatten(), bins=u_bins)
            a, _, _ = scipy.stats.binned_statistic(u_flat, F_D_sq_this[t_index, :, :].flatten(), statistic='count', bins=u_bins)
        
            # it would be quicker for F_D_sq to have dimensions u_x, u_y instead of |u|, and then we could do the binned_statistic
            # afterwards, but this would make the size of the array huge

        # we calculate these now, for the callback (otherwise we could calc them after the loop finishes)
        F_D_sq_avg = np.nanmean(F_D_sq, axis=0) # average over time origins
        num_points = np.count_nonzero(~np.isnan(F_D_sq), axis=0) # count num data points at each time origin
        F_D_sq_std = np.nanstd (F_D_sq, axis=0) # average over time origins
        F_D_sq_unc = F_D_sq_std / num_points

        u_avg = ( u_bins[:-1] + u_bins[1:] ) / 2
        
        k_avg = 2 * np.pi * u_avg
        k_x = 2 * np.pi * u_x
        k_y = 2 * np.pi * u_y

        callback(time_origin_index, k_avg, F_D_sq_avg, F_D_sq_unc, used_times*time_step, F_D_sq, time_origins, F_D_sq_noradial_sum/noradial_sum_i, k_x, k_y)

    print('final nanfrac', common.nanfrac(F_D_sq_avg))
        
    return k_avg, used_times*time_step, F_D_sq_avg, F_D_sq_unc, use_every_nth_frame, F_D_sq, time_origins, F_D_sq_noradial_sum/noradial_sum_i, k_x, k_y