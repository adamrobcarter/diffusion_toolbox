import intensity_correlation.intensity_correlation as intensity_correlation
import numpy as np
import matplotlib.pyplot as plt
import time
import common

for file in common.files_from_argv('preprocessing/data', 'stack_'):

    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    pixel_size = data['pixel_size']

    downsample = 4
    stack = stack[:, ::downsample, ::downsample]

    t0 = time.time()

    I_corr, I_mean, bins = intensity_correlation.go(stack, pixel_size)
    common.save_data(f'intensity_correlation/data/correlated_{file}.npz', 
                     I_corr=I_corr, I_mean=I_mean, bins=bins)

    # np.savez(f'data/correlated_{file}.npz', box_sizes=box_sizes, counted_intensity_diffs=counted_intensity_diffs, 
    #          avg_intensities=avg_intensities, pixel_size=pixel_size, computation_time=time.time()-t0)