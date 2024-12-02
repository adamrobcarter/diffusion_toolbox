import countoscope_old as countoscope
import numpy as np
import time
import sys
import common


if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles     = data['particles']
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']
        
        window_size = min(common.get_used_window(file, window_size_x, window_size_y))
        if '066' in file:
            num_boxes = 100
        else:
            num_boxes = 30

        if file.startswith('eleanorlong'):
            pixel_size    = data['pixel_size']
            box_sizes = np.logspace(np.log10(0.288/2), np.log10(0.9*288), num_boxes) # N was 35, but 70 for eleanorlong066
            # print('aaaa', 0.8*window_size/pixel_size)
            # box_sizes_px = np.array([0.9*window_size/pixel_size])
            # sep_sizes = 17 - box_sizes
            # sep_sizes = 9 - box_sizes # moreoverlap

        elif file.startswith('marine'):
            box_sizes = np.logspace(np.log10(0.2), np.log10(0.9*window_size), 10)
            sep_sizes = 17 - box_sizes
        elif file.startswith('marine2'):
            box_sizes_px = np.array([1,  2,  4,  8,  16,  32])
            sep_sizes_px = 7 - box_sizes_px

        elif file.startswith('sim_') or file.startswith('brennan'):
            box_sizes = np.logspace(np.log10(0.288/2), np.log10(0.9*288), num_boxes)
            sep_sizes = 7 - box_sizes # moremoreoverlap

        else:
            box_sizes = np.array([1, 2, 4, 8])
            sep_sizes = 100-box_sizes

        print('largest box', box_sizes[-1])


        output_filename = f'box_counting/data/counted_{file}'
        # output_filename += '_moreoverlap'
        output_filename += '.npz'

        t0 = time.time()
        
        for spacings in [7, 10, 14, 19, 25, 31, 38]:
            sep_sizes = 7 - box_sizes # moremoreoverlap
            output_filename = f'box_counting/data/counted_{file}'
            calc_and_save(box_sizes, sep_sizes, data, particles,
                output_filename, save_counts=False,
                save_data=True)
        print(f'took {time.time()-t0:.0f}s')