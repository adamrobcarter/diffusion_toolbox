# import intensity_counting.intensity_countoscope as intensity_countoscope
import sys
sys.path.append('/home/acarter/intensity_countoscope')
import intensity_countoscope
import numpy as np
import time
import common
import tqdm

print(intensity_countoscope.__file__)

files = ['sim', 'sim_downsampled', 'alice_0.02', 'alice_0.34', 'alice_0.66']
files = ['eleanor', 'sim', 'alice_0.02']

EVERY = 2

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    time_step         = data['time_step']
    particle_diameter = data.get('particle_diameter')

    assert np.all(stack >= 0)

    progress = tqdm.tqdm(total=stack.shape[0]*(stack.shape[0]-1)//2, desc='checking for identical frames')
    for i in range(stack.shape[0]):
        for j in range(i+1, stack.shape[0]):
            diff = np.abs(stack[i, :, :] - stack[j, :, :])
            assert diff.sum() > 0
            if diff.sum() < 100000:
                print(diff.sum())
            # assert not np.allclose(stack[i, :, :], stack[j, :, :]), f'frames {i} and {j} are identical'
            progress.update()
    # if file == 'eleanorlong':
    #     pixel_size = 0.17 # so the boxes are the same size as aliceXXX

    t0 = time.time()


    box_sizes_px = np.array([ 1,  2,  4,  8, 16, 32])
    sep_sizes_px = 60 - box_sizes_px

    box_sizes = box_sizes_px * pixel_size
    sep_sizes = sep_sizes_px * pixel_size

    # print('subtracting mean')
    # stack = stack - stack.mean(axis=0)

    stack = stack[::EVERY, :, :]

    box_sizes, counted_intensity_diffs, avg_intensities, variances, all_counts = intensity_countoscope.go(stack, box_sizes, sep_sizes, pixel_size)

    filename = f'intensity_counting/data/counted_{file}'
    if EVERY != 1:
        filename += f'_every{EVERY}'
    filename += '.npz'
    common.save_data(filename,
                     box_sizes=box_sizes, sep_sizes=sep_sizes, counted_intensity_diffs=counted_intensity_diffs, 
                     avg_intensities=avg_intensities, variances=variances, pixel_size=pixel_size, time_step=time_step,
                     particle_diameter=particle_diameter, computation_time=time.time()-t0)
#     np.savez(f'intensity_counting/data/all_counts_{file}.npz', box_sizes=box_sizes, counts=all_counts, 
#             avg_intensities=avg_intensities, variances=variances, pixel_size=pixel_size, time_step=time_step,
#             particle_diameter=particle_diameter, computation_time=time.time()-t0)