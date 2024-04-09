import numpy as np
import common
import scipy.stats
import tqdm
import numba

import time

@numba.njit
def do_binning(slice_length, F_D_sq_this, u_flat, u_bins):
    F_D_sq = np.zeros((slice_length, len(u_bins)-1))
    for t_index in numba.prange(slice_length):
        # assert u.shape == F_D_sq_this[t_index, :, :].shape
        F_D_sq[t_index, :], _, _ = common.numba_binned_statistic(u_flat, F_D_sq_this[t_index, :, :].flatten(), bins=u_bins)
    return F_D_sq

max_K = 10

def calc(stack, pixel_size, time_step, num_k_bins):
    I = stack # we use Cerbino's notation

    use_every_nth_frame = max(int(stack.size / 1e8), 1)
    print(f'automatic: use every {use_every_nth_frame}th frame')
    # use_every_nth_frame = 1 # 10 is a good number for Alice, 1 for Marine]

    time_origins = range(0, I.shape[0], use_every_nth_frame)

    used_times = common.exponential_integers(1, I.shape[0]-1) - 1

    min_K = 2*np.pi/( min(stack.shape[1], stack.shape[2]) * pixel_size )
    k_bins = np.logspace(np.log10(min_K), np.log10(max_K), num_k_bins)
    u_bins = k_bins / (2*np.pi)
    num_u_bins = u_bins.size - 1 # minus 1 because we specify right and left of final bin
    
    F_D_sq = np.full((len(time_origins), used_times.size, num_u_bins), np.nan)
    print('F_D_sq array size', common.arraysize(F_D_sq))

    # a = 0
    # b = 0

    for time_origin_index, time_origin in tqdm.tqdm(enumerate(time_origins), total=len(time_origins)):
        
        # t0 = time.time()

        time_indexes = time_origin+used_times
        time_indexes = time_indexes[time_indexes < I.shape[0]] # remove any that are longer than the data
        slice_length = time_indexes.size

        D = I[time_indexes, :, :] - I[time_origin, :, :] # need -1 cause used_times[0]==1

        u_x, u_y, F_D = common.fourier_2D(D, spacing=pixel_size, axes=(1, 2))
        
        F_D_sq_this = np.abs(F_D)**2
        # print('how real:', np.nanmean(np.imag(F_D))/np.nanmean(np.real(F_D)), np.nanmean(np.imag(F_D)), np.nanmean(np.real(F_D)))
    
        u = np.sqrt(u_x**2 + u_y**2)

        # t1 = time.time()

        # this loop below is very surely the time bottleneck
        u_flat = u.flatten()
        for t_index in range(slice_length):
            F_D_sq[time_origin_index, t_index, :], _, _ = scipy.stats.binned_statistic(u_flat, F_D_sq_this[t_index, :, :].flatten(), bins=u_bins)
            # should we check that this line is doing exactly what we expect?
        # F_D_sq[time_origin_index, :, :] = do_binning(slice_length, F_D_sq_this, u_flat, u_bins)
            
        # t2 = time.time()

        # a += t1 - t0
        # b += t2 - t1

    # print(f'{a:.3f}:{b:.3f}')

    print('max u', u.max(), u_bins.max())

    F_D_sq_avg = np.nanmean(F_D_sq, axis=0) # average over time origins

    # assert np.isnan(F_D_sq_avg).sum()/F_D_sq_avg.size < 0.1, f'F_D_sq_avg was {np.isnan(F_D_sq_avg).sum()/F_D_sq_avg.size:.2f} nan'

    u_avg = ( u_bins[:-1] + u_bins[1:] ) / 2
    
    k_avg = 2 * np.pi * u_avg

    return k_avg, used_times*time_step, F_D_sq_avg, use_every_nth_frame