import countoscope.Numba_Box_Count_Stats as countoscope
import numpy as np
import time
import sys
import common

t0 = time.time()

for file in sys.argv[1:]:
    a = 1.395 #radius of particles

    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles         = data['particles']
    time_step         = data['time_step']
    pixel_size        = data['pixel_size']
    num_timesteps     = data['num_timesteps']
    particle_diameter = data['particle_diameter']
    # particles[:, [2]] += 1 # make t 1-based
    Nframes = None
    window_size_x = None
    window_size_y = None

    if file == 'eleanorlong':
        pixel_size = 0.17 # so the boxes are the same size as aliceXXX
    
    box_sizes = np.array([1,   2,   4,   8,   16,  32,  64,  128]) * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    sep_sizes = np.array([4*a, 4*a, 3*a, 2*a, 2*a, 2*a, 2*a, 2*a])
    # bigger sep at low box sizes to reduce the number of boxes, as we already have loads of stats for small boxes
    # sep[0] = -60

    N2_mean, N2_std, N_stats = countoscope.Calc_and_Output_Stats(data='raw_data/0.02_newtrack.dat',
                                                                 window_size_x=window_size_x, window_size_y=window_size_y, 
                                                                 box_sizes=box_sizes,
                                                                #  box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
                                                                 sep_sizes=sep_sizes,)
    

    # also let's calc the packing fraction now as it'll be useful in the future
    window_width  = particles[:, 0].max() - particles[:, 0].min()
    window_height = particles[:, 1].max() - particles[:, 1].min()
    density = particles.shape[0]/num_timesteps / (window_width * window_height)
    pack_frac = np.pi/4 * density * particle_diameter**2

    # np.savez(f'box_counting/data/counted_{file}.npz', N2_mean=N2_mean, N2_std=N2_std,
    #          N_stats=N_stats, box_sizes=box_sizes, sep_sizes=sep_sizes,
    #          time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter)
    
t1 = time.time()
print(f'done in {t1-t0:.0f}s')