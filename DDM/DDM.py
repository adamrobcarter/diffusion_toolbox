import numpy as np
import common
import scipy.stats
import tqdm

def calc(stack, pixel_size, time_step):
    I = stack # we use Cerbino's notation
    print('I', I.shape)

    use_every_nth_frame = 10
    time_origins = range(0, I.shape[0], use_every_nth_frame)
    # time_origins = [0]

    print((len(time_origins), *I.shape))

    used_times = common.exponential_integers(1, 1023)
    
    F_D_sq = np.full((len(time_origins), used_times.size, I.shape[1], I.shape[2]), np.nan)
    print(common.arraysize(F_D_sq))

    for time_origin_index, time_origin in tqdm.tqdm(enumerate(time_origins), total=len(time_origins)):
        time_indexes = time_origin+used_times-1
        time_indexes = time_indexes[time_indexes < I.shape[0]] # remove any that are longer than the data
        slice_length = time_indexes.size

        D = I[time_indexes, :, :] - I[time_origin, :, :] # need -1 cause used_times[0]==1

        u_x, u_y, F_D = common.fourier_2D(D, spacing=pixel_size, axes=(1, 2)) # is this F_D or F_D^2?
        
        F_D_sq[time_origin_index, :slice_length, :, :] = np.abs(F_D)**2
        print(u_x.shape, u_y.shape, F_D_sq.shape)
        u = np.sqrt(u_x**2 + u_y**2)

    F_D_sq = np.nanmean(F_D_sq, axis=0) # average over time origins
    # print(u_x[0], u[0], u_x[1], u[1])
    num_bins = 250
    
    u_bins = np.linspace(0, u.max(), num_bins+1)
    F_D_avg = np.full((F_D_sq.shape[0], num_bins), np.nan)

    for t in range(F_D_sq.shape[0]):
        assert u.shape == F_D_sq[t, :, :].shape
        F_D_avg[t, :], _, _ = scipy.stats.binned_statistic(u.flatten(), F_D_sq[t, :, :].flatten(), bins=u_bins)

    u_avg = ( u_bins[:-1] + u_bins[1:] ) / 2
    
    k_avg = 2 * np.pi * u_avg

    return k_avg, used_times*time_step, F_D_avg