import numpy as np
import numba
import common
import tqdm


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


def reshape(particles):    
    # reformat to be (num particles) x (num timesteps) x (x, y)
    
    num_particles = int(particles[:, 3].max()) + 1
    num_timesteps = int(particles[:, 2].max()) + 1

    assert particles[:, 3].min() == 0
    assert particles[:, 2].min() == 0

    data = np.full((num_particles, num_timesteps, 2), np.nan)

    for row_index in tqdm.trange(particles.shape[0]):
        x, y, t, num = particles[row_index]
        data[int(num)-0, int(t)-0, :] = x, y # -0 cause they are 0-based

    return data

def calc_internal(data):
    x_MSDs = msd_matrix(data[:, :, 0])
    y_MSDs = msd_matrix(data[:, :, 1])
    MSDs = x_MSDs + y_MSDs

    mean, std = np.nanmean(MSDs, axis=0), np.nanstd(MSDs, axis=0)

    assert common.nanfrac(mean) < 0.1, f'MSDs nanfrac was {common.nanfrac(mean)}'

    return mean, std

def calc(particles):
    data = reshape(particles)
    # data should be of shape (num particles) x (num timesteps) x (2)

    return calc_internal(data)

def calc_individuals(particles):
    data = reshape(particles)
    # data should be of shape (num particles) x (num timesteps) x (2)
    
    x_MSDs = msd_matrix(data[:, :, 0])
    y_MSDs = msd_matrix(data[:, :, 1])
    MSDs = x_MSDs + y_MSDs

    return MSDs

def calc_centre_of_mass(data, groupsize):

    np.random.shuffle(data) # shuffles along first axes (num particles)

    groups = np.array_split(data, np.ceil(data.shape[0]/groupsize), axis=0) # ceil needed otherwise the groups will be 1 too big because of roundind down
    
    centre_of_masses = np.array([np.nanmean(group, axis=0) for group in groups])
    print('need to consider that many of your group might be nan')
    return calc_internal(centre_of_masses)


"""
@numba.njit()
def calc_msd_for_shift(data, shift_size, shift_length, exponent, shift_index):
    # data should be of shape (num particles) x (num timesteps) x (2)
    assert data.shape[0] > 1 # num_particles == 1 fails because the dimension gets squeezed
    # t0 = time.time()
    
    shift = shift_index * shift_size
    assert shift+shift_length <= data.shape[1], 'shift_index * shift_size + shift_length was longer than num_timesteps'
    
    data_to_consider = data[:, shift:shift+shift_length, :]
    # below loop is numba compatible replacement of this line:
    # is_nan = np.any(np.isnan(data_to_consider), axis=(1, 2))
    is_nan = np.full((data_to_consider.shape[0]), True)
    for particle_i in range(0, data_to_consider.shape[0]):
        is_nan[particle_i] = np.any(np.isnan(data_to_consider[particle_i, :, 0])) # only test x-value

    #print(f"dropping {is_not_nan.sum()/data_to_consider.shape[0]*100:.0f}%")
                        
    data_to_use = data_to_consider[~is_nan, :, :]

    r0 = data_to_use[:, 0, :]
    assert r0.size != 0
    assert np.isnan(r0).sum() == 0
    r0 = r0[:, np.newaxis, :] # could remove this line if r0 = data_to_use[:, [0], :]

    # with viztracer.get_tracer().log_event("test1"):
    displacements = (data_to_use - r0) # r0 shape is broadcasted up to data_to_use shape

    # with viztracer.get_tracer().log_event("test2"):
    square_displacements = displacements[:, :, 0]**2 + displacements[:, :, 1]**2

    # square_displacements = square_displacements ** (exponent/2)

    # below loop is numba-compatibe equivalent of this line:
    # msd = square_displacements.mean(axis=0)
    msd = np.full((shift_length), np.nan)
    for time_i in range(0, shift_length):
        msd[time_i] = 1
        msd[time_i] = square_displacements[:, time_i].mean()

    #mean_displacement = displacements.mean(axis=0)
    #print(mean_displacement.shape)
    #mean_displacement_length = np.sqrt(mean_displacement[:, 0]**2 + mean_displacement[:, 1]**2)

    # msds[shift_index, :] = msd
    #mds [shift_index, :] = mean_displacement_length
    # print('done')
    # progress.tick()
    # t1 = time.time()
    # print(f'{(t1 - t0)*100}ms')
    return msd

# @viztracer.log_sparse(stack_depth=3)
# @viztracer.log_sparse()
@numba.njit(parallel = True)
# @numba.njit()
def MSD(data, num_shifts, shift_length, time_step, exponent=2):
    # data should be of shape (num particles) x (num timesteps) x (2)
    num_timesteps = data.shape[1]

    assert num_shifts > 0
    assert shift_length > 0

    shift_size = int((num_timesteps - shift_length) / num_shifts)
    
    assert shift_size > 0
    shift_overlap = (shift_length-shift_size)/shift_length
    shift_overlap = shift_overlap if shift_overlap > 0 else 0
    # print(f"{num_shifts} shifts of size {shift_size} (length = {shift_length}, overlap = {shift_overlap*100:.2f}%)")
    print(num_shifts, 'shifts of size', shift_size, '(length =', shift_length, 'overlap =', shift_overlap, ')')

    msds = np.full((num_shifts, shift_length), np.nan)# * units.micrometer**exponent

    # we could get extra data at low t if we wanted, at the moment we average over the same number
    # regardless of how big or small t is
    
    parallel = False

    print('starting main MSD section')

    if not parallel:
        # progress = Progress(num_shifts)
        # times = []
        for shift_index in numba.prange(0, num_shifts):
            c = calc_msd_for_shift(data, shift_size, shift_length, exponent, shift_index)
            msds[shift_index, :] = c# * units.micrometer**2
            # times.append(c[1])
            # progress.tick()
            print(shift_index, '/', num_shifts)
        # print(f'time: {np.mean(times)*1000:.0f} +- {np.std(times)*1000:.0f} ms')
        # progress.done()

    # else:
    #     chunksize = 9
    #     print(f'{multiprocessing.cpu_count()} cores, {math.ceil(num_shifts/chunksize)} chunks of size {chunksize}')
    #     print(f'data size: {arraysize(data)}')
    #     with multiprocessing.Pool() as pool:

    #         # def task_done(result):
    #         #     progress.tick()
    #         #     print('tick', flush=True)

    #         # shared_array_base = multiprocessing.Array(ctypes.c_double, int(np.prod(data.shape)))
    #         # data_shared = np.ctypeslib.as_array(shared_array_base.get_obj())
    #         # data_shared = data_shared.reshape(*data.shape)
    #         # data_shared[:, :, :] = data

    #         task = functools.partial(calc_msd_for_shift, data, shift_size, shift_length, exponent)

    #         list_of_msds = tqdm.tqdm(pool.imap(task, range(0, num_shifts), chunksize=chunksize), total=num_shifts)
    #         #print(type(list_of_msds))

    #         list_of_msds = list(list_of_msds)
    #         # print(type(a))
    #         # print(type(a[0]), a[0])
    #         # b = np.array(list_of_msds)
    #         # print(type(b), b)
    #         # print(type(list_of_msds))
    #         # print(type(list_of_msds[0]))

    #         # while not task.ready()

    #         # list_of_msds = tasks.get()
    #         times = []

    #         # now turn that list of np arrays back into a 2d np array
    #         for i in range(0, num_shifts):
    #             msds[i, :] = list_of_msds[i][0] * units.micrometer**2
    #             times.append(list_of_msds[i][1])

    #         print(f'time: {np.mean(times)*1000:.0f} +- {np.std(times)*1000:.0f} ms')
            


    #mean_msd = msds.mean(axis=0)
    #mean_md  = mds .mean(axis=0)
            
    print('finished main MSD section')

    mean_msd = np.full((shift_length), np.nan)#~np.any(np.isnan(data_to_consider), axis=(1, 2))  # !viztracer: log
    for time_i in numba.prange(0, shift_length):
        mean_msd[time_i] = msds[:, time_i].mean()

    print('finished means')

    t = np.arange(0, mean_msd.size) * time_step

    t        = t       [~np.isnan(mean_msd)]
    mean_msd = mean_msd[~np.isnan(mean_msd)]

    return t, mean_msd
    """