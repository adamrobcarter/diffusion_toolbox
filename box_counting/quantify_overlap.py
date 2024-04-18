import common
import box_counting.count
import numpy as np

NUM_SPLITS = 10

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    # box_sizes_px = np.array([ 1,  2,  4,  8, 16, 32, 64, 128])
    # # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    # sep_sizes_px = np.array([50, 50, 50, 50, 50, 50, 50,  50])
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    max_t = particles[:, 2].max()

    for split_i in range(NUM_SPLITS):
        t_low  =  split_i      * (max_t // NUM_SPLITS)
        t_high = (split_i + 1) * (max_t // NUM_SPLITS)

        t_indexes = (t_low <= particles[:, 2]) & (particles[:, 2] < t_high)
        particles_t = particles[t_indexes, :]

        box_sizes_px = np.array([1,  16,  256,])
        # box_sizes_px = np.array([25])

        for spacing in [100.5, 30.5, 10.5, 3.5]:
            output_filename = f'box_counting/data/counted_{file}_qo_split{split_i}_spacing{spacing}.npz'
            sep_sizes_px = spacing - box_sizes_px
            box_counting.count.calc_and_save(box_sizes_px, sep_sizes_px, data, particles_t, output_filename)