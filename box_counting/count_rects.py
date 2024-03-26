import countoscope
import numpy as np
import time
import sys
import common

def calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, output_file_name, drift):
    t0 = time.time()

    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles                  = data['particles']
    time_step                  = data['time_step']
    pixel_size                 = data['pixel_size']
    num_timesteps              = data['num_timesteps']
    particle_diameter          = data['particle_diameter']
    # particle_diameter_calced   = data['particle_diameter_calced']
    particle_diameter_calced   = data.get('particle_diameter_calced')
    # particles[:, [2]] += 1 # make t 1-based
    window_size_x = None
    window_size_y = None

    if file == 'eleanorlong':
        pixel_size = 0.17 # so the boxes are the same size as aliceXXX

    added_drift_x = 0
    added_drift_y = 0
    if drift:
        added_drift_x = 0.03

        particles = common.add_drift(particles, added_drift_x, added_drift_y)
    
    # box_sizes_x = box_sizes_x_px * pixel_size 
    # box_sizes_y = box_sizes_y_px * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    # sep_sizes = sep_sizes_px * pixel_size 
    # bigger sep at low box sizes to reduce the number of boxes, as we already have loads of stats for small boxes
    # sep[0] = -60

    N2_mean, N2_std, N_stats, counts = countoscope.calculate_nmsd(data=particles,
                                                                 window_size_x=window_size_x, window_size_y=window_size_y,
                                                                 box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
                                                                 sep_sizes=sep_sizes,)
    

    # also let's calc the packing fraction now as it'll be useful in the future
    window_width  = particles[:, 0].max() - particles[:, 0].min()
    window_height = particles[:, 1].max() - particles[:, 1].min()
    density = particles.shape[0]/num_timesteps / (window_width * window_height)
    pack_frac = np.pi/4 * density * particle_diameter**2

    np.savez(output_file_name, N2_mean=N2_mean, N2_std=N2_std,
             N_stats=N_stats, box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y, sep_sizes=sep_sizes,
             time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
             particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0,
             added_drift_x=added_drift_x/time_step, added_drift_y=added_drift_y/time_step)
    # np.savez(f'box_counting/data/all_counts_{file}.npz', counts=counts,
    #          N_stats=N_stats, box_sizes=box_sizes, sep_sizes=sep_sizes,
    #          time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
    #          particle_diameter_calced=particle_diameter_calced)

if __name__ == '__main__':
    for file in sys.argv[1:]:
        box_sizes = np.array([0.4,  0.8, 1.6, 3.2, 6.4, 12.8, 25.6])
        # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
        sep_sizes = np.array([8.0, 7.6, 6.8, 5.2, 2, 2, 2])

        box_sizes_x = box_sizes
        box_sizes_y = box_sizes
        filename = f'box_counting/data/counted_{file}_rects_nodrift.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, False)

        filename = f'box_counting/data/counted_{file}_rects_eq.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, True)

        box_sizes_x = box_sizes
        box_sizes_y = 130
        filename = f'box_counting/data/counted_{file}_rects_tall.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, True)

        box_sizes_x = 130
        box_sizes_y = box_sizes
        filename = f'box_counting/data/counted_{file}_rects_wide.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, True)

        filename = f'box_counting/data/counted_{file}_rects_wide_nodrift.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, False)

        box_sizes_x = np.array([1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4])
        box_sizes_y = box_sizes_x / 8
        filename = f'box_counting/data/counted_{file}_rects_wide_aspect.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, False)

        filename = f'box_counting/data/counted_{file}_rects_wide_aspect_drift.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, True)

        box_sizes_y = np.array([1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4])
        box_sizes_x = box_sizes_y / 8
        filename = f'box_counting/data/counted_{file}_rects_tall_aspect.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, False)

        filename = f'box_counting/data/counted_{file}_rects_tall_aspect_drift.npz'
        calc_and_save(box_sizes_x, box_sizes_y, sep_sizes, file, filename, True)