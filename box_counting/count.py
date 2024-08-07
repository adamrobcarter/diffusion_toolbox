import countoscope
import numpy as np
import time
import sys
import common

def calc_and_save(box_sizes_px, sep_sizes_px, data, particles, output_file_name, extra_to_save={}, save_counts=False, use_old_overlap=False):
    t0 = time.time()

    a = 1.395 #radius of particles

    # data = common.load(f'particle_detection/data/particles_{file}.npz')
    # particles                = data['particles']
    time_step                = data['time_step']
    pixel_size               = data.get('pixel_size', 1)
    num_timesteps            = data['num_timesteps']
    particle_diameter        = data.get('particle_diameter', np.nan)
    particle_diameter_calced = data.get('particle_diameter_calced')
    depth_of_field           = data.get('depth_of_field')
    window_size_x            = data.get('window_size_x')
    window_size_y            = data.get('window_size_y')

    # if file == 'eleanorlong':
    #     pixel_size = 0.17 # so the boxes are the same size as aliceXXX
    
    box_sizes = box_sizes_px * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    sep_sizes = sep_sizes_px * pixel_size 
    # bigger sep at low box sizes to reduce the number of boxes, as we already have loads of stats for small boxes
    # sep[0] = -60

    results = countoscope.calculate_nmsd(data=particles,
                                         window_size_x=window_size_x, window_size_y=window_size_y, 
                                         box_sizes=box_sizes,
                                         #  box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
                                         sep_sizes=sep_sizes,
                                         use_old_overlap=use_old_overlap)
    
    N2_mean = results.nmsd
    N2_std = results.nmsd_std
    # N2_std = results.n, N_stats, counts

    if save_counts:
        extra_to_save['counts'] = results.counts

    # print(f'drift: {common.find_drift(particles)}')
 

    # also let's calc the packing fraction now as it'll be useful in the future
    window_width  = particles[:, 0].max() - particles[:, 0].min()
    window_height = particles[:, 1].max() - particles[:, 1].min()
    density = particles.shape[0]/num_timesteps / (window_width * window_height)
    pack_frac = np.pi/4 * density * particle_diameter**2

    common.save_data(output_file_name, N2_mean=N2_mean, N2_std=N2_std,
             box_sizes=box_sizes, sep_sizes=sep_sizes,
             num_boxes=results.num_boxes, N_mean=results.N_mean, N_var=results.N_var,
             N_var_mod=results.N_var_mod, N_var_mod_std=results.N_var_mod_std, N_mean_std=results.N_mean_std,
             time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
             particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0,
             depth_of_field=depth_of_field, old_overlap=use_old_overlap,
             box_coords=results.box_coords,
             window_size_x=window_size_x, window_size_y=window_size_y,
             **extra_to_save)
    # np.savez(f'box_counting/data/raw_counts_{file}.npz', counts=counts,
    #          N_stats=N_stats, box_sizes=box_sizes, sep_sizes=sep_sizes,
    #          time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
    #          particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0,
    #          depth_of_field=depth_of_field, **extra_to_save)

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        # box_sizes_px = np.array([ 1,  2,  4,  8, 16, 32, 64, 128])
        # # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
        # sep_sizes_px = np.array([50, 50, 50, 50, 50, 50, 50,  50])

        # box_sizes_px = np.array([ 2,  4,  8, 16, 32])
        # sep_sizes_px = np.array([20, 20, 20, 20, 20])
        # box_sizes_px = np.array([128, 128])
        # sep_sizes_px = np.array([10,  -10])
        # sep_sizes_px = np.array([10,  -10])

        # for doing timescale integral:
        box_sizes_px = np.array([1,  2,  4,  8,  16,  32,  64, 128, 256])
        sep_sizes_px = 50 - box_sizes_px
        
        box_sizes_px = np.array([1,  2,  4,  8,  16,  32,  64, 128, 256])
        sep_sizes_px = np.array([49,  48,  46,  42,  24,  8,  8, 8, 8])

        # box_sizes_px = np.array([128, 128, 128, 128, 128, 128, 128])
        # sep_sizes_px = np.array([ 10, -10, -30, -50, -70, -90,-110])

        # box_sizes_px = np.array([128, 128, 128, 128, 128, 128, 128])
        # sep_sizes_px = np.array([ 10, -10, -30, -50, -70, -90,-110])

        # box_sizes_px -= 60
        # sep_sizes_px += 60

        # sep_sizes_px = np.linspace(50, -50, 41) # 41
        # sep_sizes_px = np.linspace(-120, -100, 11)
        # sep_sizes_px = np.linspace(10, -110, 7) # 41
        # box_sizes_px = np.full_like(sep_sizes_px, 128)

        sep_sizes_px = np.linspace(5, -5, 2) # 23
        sep_sizes_px = [-13]
        # sep_sizes_px = [20, -20]
        box_sizes_px = np.full_like(sep_sizes_px, 32)


        if file.startswith('marine'):
            box_sizes_px = box_sizes_px[1:]
            sep_sizes_px = sep_sizes_px[1:]

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        # output_filename = f'box_counting/data/counted_{file}_sep0.npz'
        
        # for use_old_overlap in [True, False]:
        # # for use_old_overlap in [True, False]:
        #     ov_str = 'oldov' if use_old_overlap else 'newov'
        #     output_filename = f'box_counting/data/counted_{file}_plateaus_{ov_str}.npz'
        #     calc_and_save(box_sizes_px, sep_sizes_px, data, particles,
        #         output_filename, save_counts=False, use_old_overlap=use_old_overlap)

        # for spacing in [9, 13, 17, 23, 31, 51]:
        # box_sizes_px = np.array([1, 2, 4, 8, 16, 32, 64])
        # sep_sizes_px = spacings_px - box_sizes_px
        # # spacing = 9
        # # sep_sizes_px = np.array([8, 7, 5, 1, -7, -23, -55])
        # spacings_px  = np.array([9, 9, 9, 9, 9, 4.5, 3.5])
        # sep_sizes_px = spacings_px - box_sizes_px
        # sep_sizes_px = np.array([8, 7, 5, 1, -7, -26, -59])

        # box_sizes_px = box_sizes_px[::-1]
        # sep_sizes_px = sep_sizes_px[::-1]
        # sep_sizes_px = np.full_like(box_sizes_px, 20)


        # for eleanorlong timescaleint
        # getting rid of 160 to make the presentation nicer
        box_sizes_px = np.logspace(np.log10(0.5), np.log10(100), 27)
        sep_sizes_px = 69 - box_sizes_px

        # for alice0.66 timescaleint        
        # box_sizes_px = np.array([1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64])

        # for strange brennan rises
        # box_sizes_px = np.logspace(np.log10(1), np.log10(160), 20)
        # sep_sizes_px = 69 - box_sizes_px
        # sep_sizes_px[sep_sizes_px < 2] = 2

        output_filename = f'box_counting/data/counted_{file}.npz'


        t0 = time.time()
        calc_and_save(box_sizes_px, sep_sizes_px, data, particles,
            output_filename, save_counts=False, use_old_overlap=False)
        print(f'took {time.time()-t0:.0f}s')