import os
import tifffile
import numpy as np
import datetime
import scipy.fft

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