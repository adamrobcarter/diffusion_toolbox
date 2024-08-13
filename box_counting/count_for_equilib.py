import time
import common
from box_counting.count import calc_and_save

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        
        box_sizes_px = [200]
        sep_sizes_px = [0]

        output_filename = f'box_counting/data/counted_{file}_for_equilib.npz'

        t0 = time.time()
        calc_and_save(box_sizes_px, sep_sizes_px, data, particles,
            output_filename, save_counts=True, use_old_overlap=False)
        print(f'took {time.time()-t0:.0f}s')