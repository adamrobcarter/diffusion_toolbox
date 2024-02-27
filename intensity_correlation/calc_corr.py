import intensity_correlation
import numpy as np
import matplotlib.pyplot as plt
import time
import common

# data = np.load('data/stack_sim_downsampled.npz')
# stack = data['stack']

# file = 'exp'
# file = 'sim_downsampled'
# file = 'sim'
# file = 'alice_0.02'
file = 'alice_0.66'

data = common.load(f'data/stack_{file}.npz')
stack = data['stack']
pixel_size = data['pixel_size']

t0 = time.time()

I_corr, I_mean = intensity_correlation.go(stack, pixel_size)
np.savez('data/correlated_{file}.npz', I_corr=I_corr, I_mean=I_mean)

# np.savez(f'data/correlated_{file}.npz', box_sizes=box_sizes, counted_intensity_diffs=counted_intensity_diffs, 
#          avg_intensities=avg_intensities, pixel_size=pixel_size, computation_time=time.time()-t0)