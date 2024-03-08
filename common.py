import os
import tifffile
import numpy as np
import datetime
import scipy.fft
import tqdm
import matplotlib.animation
import numba

def get_directory_files(directory, extension):
    filenames = []
    all_filenames = os.listdir(directory)
    all_filenames.sort()
    for filename in all_filenames:
        if filename.endswith(extension):
            filenames.append(f'{directory}/{filename}')
    print(f'found {len(filenames)} .{extension}s in {directory}')
    return filenames

def load_tif(filename):
    tif = tifffile.imread(filename)
    return tif

def intensity_correlation(data1, data2):
    correlation = np.multiply.outer(data1, data2, dtype='int16')
    # correlation[a][b][c][d] is data1[a][b] * data2[c][d]
    print(correlation.shape, correlation.dtype)
    x = np.arange(0, data1.shape[0], dtype='float16')
    y = np.arange(0, data1.shape[1], dtype='float16')

    # x1, y1 = np.meshgrid(x, y)
    # x2, y2 = np.copy(x1), np.copy(y1)
    x1 = x[:, np.newaxis, np.newaxis, np.newaxis]
    y1 = y[np.newaxis, :, np.newaxis, np.newaxis]
    x2 = x[np.newaxis, np.newaxis, :, np.newaxis]
    y2 = y[np.newaxis, np.newaxis, np.newaxis, :]

    r_sq = (x1 - x2)**2 + (y1 - y2)**2

    print(r_sq.shape)

def load(filename):
    modified = datetime.datetime.fromtimestamp(os.path.getmtime(f'{filename}'))
    diff = datetime.datetime.now() - modified
    print(f'loading {filename}, last modified {str(diff)[:-10]} ago')

    data = np.load(f'{filename}', allow_pickle=True)

    if filename.endswith('.npz'):
        for key in data.keys():
            if data[key].shape: # array
                print(f'  loaded {key}, dtype={data[key].dtype}, shape={data[key].shape}, size={arraysize(data[key])}')
            else: # single value
                print(f'  loaded {key}, dtype={data[key].dtype}, value={data[key]}')
    else:
        print(f'  loaded, dtype={data.dtype}, shape={data.shape}, size={arraysize(data)}')
        
    return data

def arraysize(arr):
    size = arr.size * arr.itemsize
    if size < 1e3:
        # return str(size) + 'B'
        return f'{size}B'
    if size < 1e6:
        # return str(size/1e3) + 'kB'
        return f'{size/1e3:.0f}kB'
    if size < 1e9:
        # return str(size/1e6) + 'MB'
        return f'{size/1e6:.0f}MB'
    else:
        # return str(size/1e9) + 'GB'
        return f'{size/1e9:.0f}GB'

def fourier(t, x):
    print("doing fft")
    N = x.shape[0]
    r_spacing = t[1] - t[0]
    X = scipy.fft.rfft(x)
    f = scipy.fft.rfftfreq(N, r_spacing)
    # S = scipy.fft.fftshift(S)
    # k = scipy.fft.fftshift(k)
    # return k[N//2:], S[N//2:]
    X = np.abs(X)

    return f, X#f[:N//2], X[:N//2]

def r_squared(y_data, y_pred):
    # from https://stackoverflow.com/a/37899817/1103752
    residual_sum = np.sum((y_data - y_pred       )**2)
    total_sum    = np.sum((y_data - y_data.mean())**2)
    return 1 - (residual_sum / total_sum)

def save_gif(func, frames, fig, file, fps=1, dpi=300):
    progress = tqdm.tqdm(total=len(frames)+1) # unsure why +1 but it seems to be needed

    def frame(timestep):
        func(timestep)
        progress.update()

    ani = matplotlib.animation.FuncAnimation(fig, frame, frames=frames, interval=500, repeat=False)
    ani.save(file, dpi=dpi, writer=matplotlib.animation.PillowWriter(fps=fps))
    progress.close()

def famous_f(tau):
    return np.sqrt(tau / np.pi) * ( np.exp(-1/tau) - 1) + scipy.special.erf(np.sqrt(1/tau)) # countoscope eq. 2

def structure_factor_2d_hard_spheres(k, phi, sigma):
    # sigma is disk diameter
    rho = 4 * phi / (np.pi * sigma**2)

    phi = np.pi/4 * rho * sigma**2
    J0 = lambda x: scipy.special.jv(0, x)
    J1 = lambda x: scipy.special.jv(1, x)

    prefactor = np.pi / ( 6 * ( 1 - phi)**3 * k**2 )
    line1 = -5/4 * (1 - phi)**2 * k**2 * sigma**2 * J0(k * sigma / 2)**2
    line23 = 4 * ( (phi - 20) * phi + 7) + 5/4 * (1 - phi)**2 * k**2 * sigma**2
    line23factor = J1(k * sigma / 2)**2
    line4 = 2 * (phi - 13) * (1 - phi) * k * sigma * J1(k * sigma / 2) * J0(k * sigma / 2)
    c = prefactor * (line1 + line23*line23factor + line4)
    # ^^^ Thorneywork et al 2018
    
    S = 1 / (1 - rho * c) # Hansen & McDonald (3.6.10)

    return S


@numba.njit
def numba_binned_statistic(x, y, bins):
    # print(x.min(), x.max())
    # bin_edges = np.linspace(0, 100, bins+1) # +1 because you need eg. 3 points to define 2 bins
    # this will have broken if u used to supply num_bins
    bin_edges = bins
    num_bins = len(bins)-1
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

    return bin_sums / bin_counts, bin_edges, None

@numba.njit
def numba_sum_3d_axis01(array):
    # numba compatible version of array.sum(axis=(0, 1))
    sums = np.zeros((array.shape[2]))
    for i in range(array.shape[0]):
        sums[i] = array[:, :, i].sum()
    return sums

@numba.njit
def numba_nanmean_2d_axis0(array):
    nanmeans = np.zeros(array.shape[1])
    for i in range(array.shape[1]):
        nanmeans[i] = np.nanmean(array[:, i])
    return nanmeans

@numba.njit
def numba_nanstd_2d_axis0(array):
    nanstds = np.zeros(array.shape[1])
    for i in range(array.shape[1]):
        nanstds[i] = np.nanstd(array[:, i])
    return nanstds

@numba.njit
def numba_p_assert(condition, message):
    if not condition:
        print('Assertion failed: ', message)