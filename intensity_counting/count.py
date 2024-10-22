import intensity_counting.intensity_countoscope as intensity_countoscope
import numpy as np
import time
import common
import sys

files = ['sim', 'sim_downsampled', 'alice_0.02', 'alice_0.34', 'alice_0.66']
files = ['eleanor', 'sim', 'alice_0.02']

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    time_step         = data['time_step']
    particle_diameter = data.get('particle_diameter')

    # if file == 'eleanorlong':
    #     pixel_size = 0.17 # so the boxes are the same size as aliceXXX

    t0 = time.time()

    # box_sizes = (0.5, 1, 2, 4, 8, 16)
    box_sizes_px = (1, 2, 3, 4, 6, 8, 10, 16, 20, 26, 32, 64, 128)
    box_sizes_px = np.array([ 1,  2,  4,  8, 16, 32,  64, 128])
    # sep_sizes_px = (60, 50, 40, 30, 20, 20, -20, -50)
    sep_sizes_px = 60 - box_sizes_px
    # box_sizes_px = ( 1,  2,  3,  4,  6,  8, 11, 16, 22, 32, 44,  64,  84, 128)
    # sep_sizes_px = (60, 50, 50, 40, 40, 30, 30, 20, 20, 20, 14, -20, -30, -50)
    # box_sizes_px = np.unique(np.round(np.logspace(0, np.log10(128)), 50).astype('int'))
    # sep_sizes_px = np.full_like(box_sizes_px, 30)

    box_sizes = box_sizes_px * pixel_size
    sep_sizes = sep_sizes_px * pixel_size

    print('subtracting mean')
    stack = stack - stack.mean(axis=0)

    box_sizes, counted_intensity_diffs, avg_intensities, variances, all_counts = intensity_countoscope.go(stack, box_sizes, sep_sizes, pixel_size)

    common.save_data(f'intensity_counting/data/counted_{file}.npz',
                     box_sizes=box_sizes, sep_sizes=sep_sizes, counted_intensity_diffs=counted_intensity_diffs, 
                     avg_intensities=avg_intensities, variances=variances, pixel_size=pixel_size, time_step=time_step,
                     particle_diameter=particle_diameter, computation_time=time.time()-t0)
#     np.savez(f'intensity_counting/data/all_counts_{file}.npz', box_sizes=box_sizes, counts=all_counts, 
#             avg_intensities=avg_intensities, variances=variances, pixel_size=pixel_size, time_step=time_step,
#             particle_diameter=particle_diameter, computation_time=time.time()-t0)