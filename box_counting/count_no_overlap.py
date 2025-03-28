import countoscope_old as countoscope
import numpy as np
import time
import sys
import common
import box_counting.count

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']
        
        window_size = min(common.get_used_window(file, window_size_x, window_size_y))
        if '066' in file:
            num_boxes = 100
        else:
            num_boxes = 30

        if file.startswith('eleanorlong'):
            pixel_size    = data['pixel_size']
            box_sizes = np.logspace(np.log10(pixel_size/2), np.log10(0.9*window_size), num_boxes) # N was 35, but 70 for eleanorlong066


        elif file.startswith('brennan'):

            box_sizes = np.logspace(np.log10(0.288), np.log10(0.9*window_size), num_boxes)
            sep_sizes = 17 - box_sizes # moremoreoverlap

        elif file.startswith('marine'):
            box_sizes = np.logspace(np.log10(0.2), np.log10(0.9*window_size), 10)
            sep_sizes = 17 - box_sizes
        elif file.startswith('marine2'):
            box_sizes_px = np.array([1,  2,  4,  8,  16,  32])
            sep_sizes_px = 7 - box_sizes_px

        elif file.startswith('sim_'):
            box_sizes = np.logspace(np.log10(0.288), np.log10(0.9*window_size), 35)
            sep_sizes = 7 - box_sizes # moremoreoverlap

        else:
            box_sizes = np.array([1, 2, 4, 8])
            sep_sizes = 100-box_sizes

        sep_sizes = 50-box_sizes
        sep_sizes[sep_sizes < 2] = 2


        output_filename = f'box_counting/data/counted_{file}_no_overlap'
        # output_filename += '_moreoverlap'
        output_filename += '.npz'

        t0 = time.time()
        box_counting.count.calc_and_save(box_sizes=box_sizes, sep_sizes=sep_sizes, data=data,
            output_file_name=output_filename, save_counts=False,
            save_data=True)
        print(f'took {time.time()-t0:.0f}s')