import common
import box_counting.count
import numpy as np
import tqdm

NUM_SPLITS = 10
SPACINGS = [512.5, 256.5, 128.5, 64.5, 32.5, 16.5, 8.5, 4.5, 2.5, 1.5, 1, 0.5, 0.3, 0.15]
# SPACINGS = [16.5]
SPACINGS = [512.5]
SPACINGS = np.round(np.logspace(np.log10(0.4), np.log10(512.5), 30), decimals=2)
SPACINGS = np.round(np.logspace(np.log10(0.4), np.log10(512.5), 50), decimals=2)
# SPACINGS = np.round(np.logspace(np.log10(2.5), np.log10(156.5), 50), decimals=2)
print(SPACINGS)

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    # box_sizes_px = np.array([ 1,  2,  4,  8, 16, 32, 64, 128])
    # # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    # sep_sizes_px = np.array([50, 50, 50, 50, 50, 50, 50,  50])
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    max_t = particles[:, 2].max()

    # max_dim = 500
    # x_too_big = particles[:, 0] > 200
    # y_too_big = particles[:, 1] > 200
    # keep = (~x_too_big) & (~y_too_big)
    # particles = particles[keep, :]

    progress = tqdm.tqdm(total=NUM_SPLITS*len(SPACINGS))

    for split_i in range(NUM_SPLITS):
        print()
        print(f'on split {split_i+1} of {NUM_SPLITS}')
        t_low  =  split_i      * (max_t // NUM_SPLITS)
        t_high = (split_i + 1) * (max_t // NUM_SPLITS)
        print('split length', t_high-t_low)

        t_indexes = (t_low <= particles[:, 2]) & (particles[:, 2] < t_high)
        particles_t = particles[t_indexes, :]

        box_sizes_px = np.array([1, 4, 16, 64, 256])

        for spacing in SPACINGS:
            output_filename = f'box_counting/data/counted_{file}_extra4_qo_split{split_i}_spacing{spacing}.npz'
            sep_sizes_px = spacing - box_sizes_px
            box_counting.count.calc_and_save(box_sizes_px, sep_sizes_px, data, particles_t, output_filename)
            progress.update()
    
    progress.close()