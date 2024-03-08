import numpy as np
import tqdm
import numba

def autocorrFFT(x):
    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * np.conjugate(F)
    res = np.fft.ifft(PSD)
    res = (res[:N]).real  # now we have the autocorrelation in convention B
    n = N * np.ones(N) - np.arange(0, N)  # divide res(m) by (N-m)
    return res / n  # this is the autocorrelation in convention A

@numba.njit(fastmath=True)
def msd_fft1d(r):
    N = len(r)
    D = np.square(r)
    D = np.append(D, 0)
    with numba.objmode(S2='float64[:]'):
        S2 = autocorrFFT(r)
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in range(N):
        Q = Q - D[m-1] - D[N-m]
        S1[m] = Q / (N-m)
    return S1 - 2 * S2

@numba.njit(parallel=True, fastmath=True)
def msd_matrix(matrix):
    # calculates the MSDs of the rows of the provided matrix
    Nrows, Ncols = matrix.shape
    MSDs = np.zeros((Nrows,Ncols))
    for i in range(Nrows):
        #print(100.0 * ((1.0 * i) / (1.0 * N)), "percent done with MSD calc")
        MSD = msd_fft1d(matrix[i, :])
        MSDs[i,:] = MSD
    return MSDs

def go(intensities, box_sizes_um, sep_sizes_um, pixel_size):
    assert len(box_sizes_um) == len(sep_sizes_um), f'len(box_sizes_um) = {len(box_sizes_um)} != len(sep_sizes_um) = {len(sep_sizes_um)}'
    # intensity_diffs = intensities[:, :, :] - intensities[:, :, [0]]
    num_timesteps = intensities.shape[0]

    window_width  = intensities.shape[1]
    window_height = intensities.shape[2]

    # box sizes and sep come in micrometers, but we need them in pixels
    box_sizes = np.array([int(round(box_size_um / pixel_size)) for box_size_um in box_sizes_um])
    sep_sizes = np.array([int(round(sep_size_um / pixel_size)) for sep_size_um in sep_sizes_um])
    # need the round() there or sometimes going back and forth between um and pixel results in a different value

    assert np.all((box_sizes < window_width) & (box_sizes < window_height)), 'Box sizes must be smaller than window size'

    print('box sizes px', box_sizes, sep_sizes)

    msd_means, avgs, variances, all_counts = inner_loop(intensities, num_timesteps, box_sizes, sep_sizes, window_width, window_height)

    return box_sizes * pixel_size, msd_means, avgs, variances, all_counts

def inner_loop(intensities, num_timesteps, box_sizes, sep_sizes, window_width, window_height):
    msd_length = num_timesteps - 0#np.max(time_origins)
    msd_means = np.full((len(box_sizes), msd_length), np.nan)
    # msd_stds  = np.zeros((len(box_sizes), num_timesteps))
    avgs      = np.full((len(box_sizes)), np.nan)
    variances = np.full((len(box_sizes)), np.nan)

    all_counts = []

    for box_size_index in tqdm.trange(len(box_sizes)):
        box_size = box_sizes[box_size_index]
        sep      = sep_sizes[box_size_index]
        msds, avg, variance, counts = msds_for_box_size(intensities, box_size, sep, window_width, window_height, msd_length)
        msd_means[box_size_index, :] = msds.mean(axis=0)
        # msd_stds [box_size_index, :] = msd_std
        avgs     [box_size_index]    = avg
        variances[box_size_index]    = variance
        print('avg', avg)
        all_counts.append(counts)
    
    return msd_means, avgs, variances, all_counts

@numba.njit(fastmath=True, parallel=True)
def msds_for_box_size(intensities, box_size, sep, window_width, window_height, msd_length):
    # count the intensities for one box size

    num_timesteps = intensities.shape[0]

    num_boxes_x = int(np.floor(window_width  / (box_size + sep)))
    num_boxes_y = int(np.floor(window_height / (box_size + sep)))
    counts = np.full((num_timesteps, num_boxes_x * num_boxes_y), np.nan)
    print('L =', box_size, 'sep =', sep, ':', num_boxes_x, '*', num_boxes_y)
    
    box_xs = np.arange(0, num_boxes_x) * (box_size + sep) #+ sep/2
    box_ys = np.arange(0, num_boxes_y) * (box_size + sep) #+ sep/2

    if sep < 0 and box_size % np.abs(sep) == 0:
        print('Negative overlap is an exact divisor of box size. This will lead to correlated boxes.')
        

    for timestep in numba.prange(num_timesteps):
        for box_x_index in range(num_boxes_x):
            box_x = box_xs[box_x_index]
            for box_y_index in range(num_boxes_y):
                box_y = box_ys[box_y_index]

                counts[timestep, box_x_index * num_boxes_y + box_y_index] = intensities[timestep, box_x:box_x+box_size, box_y:box_y+box_size].sum()

    avg = counts.mean()
    variance = counts.var()
    # variance = np.mean(counts**2) - avg**2

    counts_t = np.transpose(counts)
    msds = msd_matrix(counts_t)
    # msd_avgs = np.mean(msds, axis=0)
    return msds, avg, variance, counts[:, 0:100:10]