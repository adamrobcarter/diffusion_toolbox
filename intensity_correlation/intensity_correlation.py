import numpy as np
import scipy.stats
import tqdm
import numba

def go(stack, pixel_size):
    return correlate_two_frames(stack[0, :, :], stack[0, :, :], pixel_size)

def correlate_two_frames(frame1, frame2, pixel_size):
    # if we fully vectorise this we exceed the available ram haha

    # xs = np.linspace(0, frame1.shape[0]*pixel_size, frame1.shape[0])

    # ys = np.linspace(0, frame1.shape[1]*pixel_size, frame1.shape[1])
    num_x_pixels = frame1.shape[0]
    num_y_pixels = frame1.shape[1]

    print(f'{num_x_pixels} * {num_y_pixels} px')


    # min_k = 2*np.pi/( min(num_x_pixels, num_y_pixels) * crop ) # note this is in pixels
    # # min_K = 2*np.pi/( min(num_x_pixels, num_y_pixels) * crop ) # note this is in pixels
    # num_k_bins = 50
    # max_k = 10 / pixel_size

    # r_bins = np.linspace(0, max_k, num_k_bins)

    # we can get a saving here because we're counting each pair of pixels in both directions
    num_bins = 50

    all_corr_binned = np.zeros((num_bins))
    # print('starting parallel computation', 'eta', 429/1.07/12)
    # for x_index_1 in tqdm.trange(0, 100):
    for x_index_1 in tqdm.trange(0, num_x_pixels):
        all_corr_binned += outer_inner(frame1, frame2, num_x_pixels, num_y_pixels, num_bins, x_index_1)
        # break
    # print('finished parallel computation')

    all_corr_binned /= (num_x_pixels * num_y_pixels)
    return all_corr_binned, np.mean([frame1.mean(), frame2.mean()])

@numba.njit(parallel=True)
def outer_inner(frame1, frame2, num_x_pixels, num_y_pixels, num_bins, x_index_1):
    all_corr_binned = np.zeros((num_bins))

    for y_index_1 in numba.prange(0, num_y_pixels):
    # for y_index_1 in numba.prange(0, 100):
        corrs, rs = inner(num_x_pixels, num_y_pixels, frame1, frame2, x_index_1, y_index_1)
        all_corr_binned += numba_binned_statistic(rs.flatten(), corrs.flatten(), num_bins)
        # break

    # corr_binned, r_binned, _ = scipy.stats.binned_statistic(rs.flatten(), corrs.flatten(), 'mean', bins=num_bins)
    
    return all_corr_binned

@numba.njit
def numba_binned_statistic(x, y, num_bins):
    # print(x.min(), x.max())
    bin_edges = np.linspace(0, 100, num_bins+1) # +1 because you need eg. 3 points to define 2 bins
    bin_sums   = np.zeros((num_bins))
    bin_counts = np.zeros((num_bins))

    for value_i in range(x.shape[0]):
        for bin_i in range(num_bins):
            # print(x[value_i], '>', bin_edges[bin_i], x[value_i] > bin_edges[bin_i])
            if x[value_i] < bin_edges[bin_i+1]:
                # print('y[value_i]', y[value_i])
                bin_sums  [bin_i] += y[value_i]
                bin_counts[bin_i] += 1
                # print(x[value_i], bin_i)
                # break

    # print('suns', bin_sums)
    # print('counts', bin_counts)
    # print()

    return bin_sums / bin_counts

@numba.njit
def inner(num_x_pixels, num_y_pixels, frame1, frame2, x_index_1, y_index_1):
    corrs = np.full((num_x_pixels, num_y_pixels), np.nan)
    r2s   = np.full((num_x_pixels, num_y_pixels), np.nan)

    xs = np.arange(0, num_x_pixels)
    ys = np.arange(0, num_y_pixels)

        # for x_index_2 in range(x_index_1, num_x_pixels):
        #     for y_index_2 in range(y_index_1, num_y_pixels):

    corrs[:, :] = frame1[x_index_1, y_index_1] * frame2[:, :]
    x2 = (x_index_1 - xs)**2
    y2 = (y_index_1 - ys)**2
    r2s  [:, :] = x2[:, np.newaxis] + y2[np.newaxis, :]

    # for x_index_2 in range(0, num_x_pixels):
    #     for y_index_2 in range(0, num_y_pixels):
    #         corrs[y_index_1, x_index_2, y_index_2] = frame1[x_index_1, y_index_1] * frame2[x_index_2, y_index_2]
    #         r2s  [y_index_1, x_index_2, y_index_2] = (x_index_1 - x_index_2)**2 + (y_index_1 - y_index_2)**2

    rs = np.sqrt(r2s)
    
    return corrs, rs