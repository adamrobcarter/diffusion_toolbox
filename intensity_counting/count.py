import intensity_counting.intensity_countoscope as intensity_countoscope
import numpy as np
import time
import common
import sys

files = ['sim', 'sim_downsampled', 'alice_0.02', 'alice_0.34', 'alice_0.66']
files = ['eleanor', 'sim', 'alice_0.02']

for file in sys.argv[1:]:

    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    time_step         = data['time_step']
    particle_diameter = data['particle_diameter']

    # if file == 'eleanorlong':
    #     pixel_size = 0.17 # so the boxes are the same size as aliceXXX

    t0 = time.time()

    # box_sizes = (0.5, 1, 2, 4, 8, 16)
    box_sizes_px = (1, 2, 3, 4, 6, 8, 10, 16, 20, 26, 32, 64, 128)
    box_sizes_px = (1, 2, 4, 8, 16, 32, 64, 128)
    box_sizes = box_sizes_px * pixel_size

    box_sizes, counted_intensity_diffs, avg_intensities = intensity_countoscope.go(stack, box_sizes, pixel_size)

    np.savez(f'intensity_counting/data/counted_{file}.npz', box_sizes=box_sizes, counted_intensity_diffs=counted_intensity_diffs, 
            avg_intensities=avg_intensities, pixel_size=pixel_size, time_step=time_step,
            particle_diameter=particle_diameter, computation_time=time.time()-t0)