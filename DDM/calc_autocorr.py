import common
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import tqdm

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']

    sequence = stack[:, 0, 0]
    autocorr = scipy.signal.correlate(sequence, sequence)

    autocorrs = np.full((autocorr.shape[0], stack.shape[1], stack.shape[2]), np.nan)

    for x in tqdm.trange(stack.shape[1]):
        for y in range(stack.shape[2]):
            sequence = stack[:, x, y] - stack[:, x, y].mean()
            # print(sequence)
            autocorrs[:, x, y] = scipy.signal.correlate(sequence, sequence)

    autocorr_avg = autocorrs.mean(axis=(1, 2))
    autocorr_std = autocorrs.std (axis=(1, 2)) / np.sqrt(stack.shape[1] * stack.shape[2])
    lags = scipy.signal.correlation_lags(sequence.shape[0], sequence.shape[0])

    common.save_data(f'DDM/data/autocorr_{file}.npz', autocorr_avg=autocorr_avg, autocorr_std=autocorr_std, t=lags*data['time_step'])