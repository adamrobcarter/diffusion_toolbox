import numpy as np
import common
import scipy.stats
import tqdm
import numba
import warnings
import time, functools
import multiprocessing

@numba.njit
def do_binning(slice_length, F_D_sq_this, u_flat, u_bins):
    F_D_sq = np.zeros((slice_length, len(u_bins)-1))
    for t_index in numba.prange(slice_length):
        # assert u.shape == F_D_sq_this[t_index, :, :].shape
        F_D_sq[t_index, :], _, _ = common.numba_binned_statistic(u_flat, F_D_sq_this[t_index, :, :].flatten(), bins=u_bins)
    return F_D_sq

max_K = 10

I = None

def binning(u_flat, u_bins, F_D_sq):
    F_D_sq_binned, _, _ = scipy.stats.binned_statistic(u_flat, F_D_sq.flatten(), bins=u_bins)
    return F_D_sq_binned 

def calc_for_time_origin(pixel_size, I, used_times, num_u_bins, u_bins, time_origin):
    time_indexes = time_origin+used_times
    time_indexes = time_indexes[time_indexes < I.shape[0]] # remove any that are longer than the data
    slice_length = time_indexes.size

    D = I[time_indexes, :, :] - I[time_origin, :, :] # need -1 cause used_times[0]==1

    u_x, u_y, F_D = common.fourier_2D(D, spacing=pixel_size, axes=(1, 2))
    # F_D = I_fourier[time_indexes, :, :] - I_fourier[time_origin, :, :]

    warnings.warn('see static_fourier.py, should you be using fftshift?')
    
    F_D_sq_this = np.abs(F_D)**2
    # print('how real:', np.nanmean(np.imag(F_D))/np.nanmean(np.real(F_D)), np.nanmean(np.imag(F_D)), np.nanmean(np.real(F_D)))

    u = np.sqrt(u_x**2 + u_y**2)

    # this loop below is very surely the time bottleneck
    u_flat = u.flatten()
    F_D_sq = np.full((used_times.size, num_u_bins), np.nan)


    # it would be quicker for F_D_sq to have dimensions u_x, u_y instead of |u|, and then we could do the binned_statistic
    # afterwards, but this would make the size of the array huge

    slices = []

    for t_index in range(slice_length):
        slices.append(F_D_sq_this[t_index, :, :])

    with multiprocessing.Pool(12) as pool:
        task = functools.partial(binning, u_flat, u_bins)
        results = pool.map(task, slices)

    for i, result in enumerate(results):
        F_D_sq[i, :] = result

    return F_D_sq

def calc(stack, pixel_size, time_step, num_k_bins):
    I = stack # we use Cerbino's notation

    # use_every_nth_frame = max(int(stack.size / 1e8), 1)
    # print(f'automatic: use every {use_every_nth_frame}th frame')
    use_every_nth_frame = 1
    warnings.warn('using every frame')
    # use_every_nth_frame = 1 # 10 is a good number for Alice, 1 for Marine]
    
    # use_every_nth_frame = max(num_timesteps / max_time_origins, 1) this is what we use in scattering_functions btw

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

    # print('doing big fourier')
    # # instead of taking the intensity differences and fourier transforming them, we could fourier transform the intensities
    # # and then take the differences. This is okay cause the fourier transform is linear. However we don't do that cause the 
    # # speed bottleneck is binned_statistic, not the fourier transform
    # u_x, u_y, I_fourier = common.fourier_2D(I, spacing=pixel_size, axes=(1, 2))
    # u = np.sqrt(u_x**2 + u_y**2)
    # print('done')

    for time_origin_index, time_origin in tqdm.tqdm(enumerate(time_origins), total=len(time_origins)):
    # with multiprocessing.Pool(12) as pool:

        task = functools.partial(calc_for_time_origin, pixel_size, I, used_times, num_u_bins, u_bins)
        result = task(time_origin)
    # results = list(tqdm.tqdm(map(task, time_origins), total=len(time_origins)))

        F_D_sq[time_origin_index, :, :] = result
    

    F_D_sq_avg = np.nanmean(F_D_sq, axis=0) # average over time origins
    F_D_sq_std = np.nanstd (F_D_sq, axis=0) # average over time origins

    # assert np.isnan(F_D_sq_avg).sum()/F_D_sq_avg.size < 0.1, f'F_D_sq_avg was {np.isnan(F_D_sq_avg).sum()/F_D_sq_avg.size:.2f} nan'

    u_avg = ( u_bins[:-1] + u_bins[1:] ) / 2
    
    k_avg = 2 * np.pi * u_avg

    return k_avg, used_times*time_step, F_D_sq_avg, F_D_sq_std, use_every_nth_frame