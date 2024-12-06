import countoscope_old as countoscope
import numpy as np
import time
import sys
import common
from box_counting.count import calc_and_save

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']
        window_size = min(window_size_x, window_size_y)
        
        window_size = min(common.get_used_window(file, window_size_x, window_size_y))

            
        box_sizes = np.logspace(np.log10(0.01*window_size), np.log10(0.9*window_size), 20)
        # print('aaaa', 0.8*window_size/pixel_size)
        # box_sizes_px = np.array([0.9*window_size/pixel_size])
        # sep_sizes = 17 - box_sizes
        # sep_sizes = 9 - box_sizes # moreoverlap
        sep_sizes = 17 - box_sizes # moremoreoverlap


        output_filename = f'box_counting/data/counted_{file}_frac_of_window.npz'

        t0 = time.time()
        calc_and_save(box_sizes=box_sizes, sep_sizes=sep_sizes, data=data,
            output_file_name=output_filename, save_counts=False,
            save_data=True)
        print(f'took {time.time()-t0:.0f}s')