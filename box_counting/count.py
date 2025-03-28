import countoscope_old as countoscope
import numpy as np
import time
import sys
import common

def calc_and_save(file, box_sizes, sep_sizes, data, output_file_name, particles=None,
                  extra_to_save={}, save_counts=False,
                  save_data=True, skip_processing=False):
    t0 = time.time()

    if particles is not None:
        assert isinstance(particles, np.ndarray)
        raise Exception('this is maybe deprecated?')

    if not save_data:
        print('WARNING: I am not going to save the data')

    a = 1.395 #radius of particles

    # data = common.load(f'particle_detection/data/particles_{file}.npz')
    # particles                = data['particles']
    time_step                = data['time_step']
    pixel_size               = data.get('pixel_size', 1)
    particle_diameter        = data.get('particle_diameter', np.nan)
    particle_diameter_calced = data.get('particle_diameter_calced')
    depth_of_field           = data.get('depth_of_field')
    window_size_x            = data.get('window_size_x')
    window_size_y            = data.get('window_size_y')
    # assert 'density' in data

    if particles == None:
        particles = data['particles']
    else:
        raise Exception('should this be deprecated, and this conditional removed?')
    print('particles', particles)

    if '_pot' in file:
        print('cropping')
        keep_x = (0 <= particles[:, 0]) & (particles[:, 0] <= window_size_x)
        keep_y = (0 <= particles[:, 1]) & (particles[:, 1] <= window_size_y)
        print(keep_x.shape, particles.shape)
        particles = particles[keep_x & keep_y, :]

    # num_timesteps = int(particles[:, 2].max()) + 1

    # if file == 'eleanorlong':
    #     pixel_size = 0.17 # so the boxes are the same size as aliceXXX
    
    # box_sizes = box_sizes_px * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    # sep_sizes = sep_sizes_px * pixel_size 
    # bigger sep at low box sizes to reduce the number of boxes, as we already have loads of stats for small boxes
    # sep[0] = -60

    # extra_params_to_countoscope = {}
    # if save_counts:
    #     extra_params_to_countoscope['return_counts'] = True

    results = countoscope.calculate_nmsd(data=particles,
                                         window_size_x=window_size_x, window_size_y=window_size_y, 
                                         box_sizes=box_sizes,
                                         #  box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
                                         sep_sizes=sep_sizes,
                                         skip_processing=skip_processing,
                                         return_counts=save_counts,
                                        #  **extra_params_to_countoscope
                                         return_histograms=False,
                                         )
    
    N2_mean = results.nmsd
    N2_std = results.nmsd_std
    # N2_std = results.n, N_stats, counts

    if save_counts:
        extra_to_save['counts'] = results.counts
        assert results.counts is not None

    # print(f'drift: {common.find_drift(particles)}')
 

    # also let's calc the packing fraction now as it'll be useful in the future
    # window_width  = window_size_x if window_size_x else particles[:, 0].max() - particles[:, 0].min()
    # window_height = window_size_y if window_size_y else particles[:, 1].max() - particles[:, 1].min()
    # density = particles.shape[0]/num_timesteps / (window_size_x * window_size_y)
    # pack_frac = np.pi/4 * density * particle_diameter**2

    output = dict(filename=output_file_name, N2_mean=N2_mean, N2_std=N2_std,
                box_sizes=box_sizes, sep_sizes=sep_sizes,
                num_boxes=results.num_boxes, N_mean=results.N_mean, N_var=results.N_var, N_var_std=results.N_var_std,
                N_var_mod=results.N_var_mod, N_var_mod_std=results.N_var_mod_std, N_mean_std=results.N_mean_std,
                N_var_time=results.N_var_time,
                N_var_sem_lb=results.N_var_sem_lb, N_var_sem_ub=results.N_var_sem_ub,
                hist_x=results.hist_x, hist_y=results.hist_y,
                time_step=time_step, particle_diameter=particle_diameter,
                particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0,
                depth_of_field=depth_of_field,
                box_coords=results.box_coords,
                pack_frac=data.get('pack_frac'), density=data.get('density', common.calc_density(particles, window_size_x, window_size_y)),
                pack_frac_given=data.get('pack_frac_given'), max_time_hours=data.get('max_time_hours'),
                window_size_x=window_size_x, window_size_y=window_size_y, pixel_size=data.get('pixel_size'),
                **extra_to_save)

    if save_data:
        common.save_data(**output)
    else:
        print('not saving data')

    return output

def get_box_spacing(file):
    if file.startswith('eleanorlong'):
        return 7

    elif file.startswith('marine'):
        return 17
    
    elif file.startswith('marine2'):
        return 7

    elif file.startswith('sim_') or file.startswith('brennan'):
        if '160' in file:
            return 7
        elif '320' in file:
            return 7
        elif '640' in file:
            return 13
        elif '1280' in file:
            return 17 # it seems this could be reduced
        else:
            raise Exception('spacing not defined')
        
    return 7

def get_box_sizes(file, data):
    window_size_x = data['window_size_x']
    window_size_y = data['window_size_y']
    
    window_size = min(common.get_used_window(file, window_size_x, window_size_y))

    if '066' in file:
        num_boxes = 100
    else:
        num_boxes = 30

    if file.startswith('eleanorlong'):
        min_box = 0.288/2
        max_box = 0.9*288

    elif file.startswith('marine'):
        raise Exception('needs reimplementing')
        box_sizes = np.logspace(np.log10(0.2), np.log10(0.9*window_size), 10)
        sep_sizes = 17 - box_sizes
    elif file.startswith('marine2'):
        raise Exception('needs reimplementing')
        box_sizes_px = np.array([1,  2,  4,  8,  16,  32])
        sep_sizes_px = 7 - box_sizes_px

    elif file.startswith('sim_') or file.startswith('brennan'):
        min_box = 0.288/2
        max_box = 0.9*288

    else:
        if 'pixel_size' in data:
            min_box = data['pixel_size']
        else:
            min_box = 1
        max_box = 0.9*window_size

    box_sizes = np.logspace(np.log10(min_box), np.log10(max_box), num_boxes)
    sep_sizes = get_box_spacing(file) - box_sizes


    print('largest box', box_sizes[-1])

    if box_sizes.max() > window_size:
        sep_sizes = sep_sizes[box_sizes < window_size]
        box_sizes = box_sizes[box_sizes < window_size]

    return box_sizes, sep_sizes

def go(file):
    data = common.load(f'particle_detection/data/particles_{file}.npz')

    box_sizes, sep_sizes = get_box_sizes(file, data)

    
    print('box sizes', box_sizes/data['particle_diameter'])
    Nboxes_nowrap = np.ceil((data['window_size_x'] - box_sizes)/(box_sizes + sep_sizes))
    Nboxes_wrap   = np.ceil((data['window_size_x']            )/(box_sizes + sep_sizes))
    print('wrapping more boxes', Nboxes_wrap/Nboxes_nowrap)

    import sys
    sys.exit()

    output_filename = f'box_counting/data/counted_{file}'
        # output_filename += '_moreoverlap'
    output_filename += '.npz'

    t0 = time.time()
    calc_and_save(file=file, box_sizes=box_sizes, sep_sizes=sep_sizes, data=data,
            output_file_name=output_filename, save_counts=False,
            save_data=True)
    print(f'took {time.time()-t0:.0f}s')

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file)