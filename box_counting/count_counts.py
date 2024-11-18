import countoscope_old as countoscope
import numpy as np
import time
import sys
import common

from box_counting.count import calc_and_save

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles     = data['particles']
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']
        
        window_size = min(window_size_x, window_size_y)


        box_sizes = np.logspace(np.log10(0.288), np.log10(0.9*window_size), 20)
        sep_sizes = 17 - box_sizes


        output_filename = f'box_counting/data/counted_counts_{file}'
        # output_filename += '_moreoverlap'
        output_filename += '.npz'

        t0 = time.time()
        calc_and_save(box_sizes, sep_sizes, data, particles,
            output_filename, save_counts=True,
            save_data=True, skip_processing=True)
        print(f'took {time.time()-t0:.0f}s')