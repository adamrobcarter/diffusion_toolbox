import countoscope.Numba_Box_Count_Stats as countoscope
import numpy as np
import time
import sys
import common

t0 = time.time()

for file in sys.argv[1:]:
    a = 1.395 #radius of particles
    # sep = np.array([2*a, 2*a, 2*a, 2*a, 2*a, 2*a, 4*a, 4*a, 4*a]) #3*a #separation between boxes
    # Box_Ls = np.array([64.0, 32.0, 16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25]) # array of box sizes to probe
    # Box_Ls = np.logspace(-2, 6, 200, base=2)

    # data = np.load(f'data/{phi}_newtrack.npz')
    # data_param = data['data']
    # window_size_x = data['window_size_x']
    # window_size_y = data['window_size_y']
    # Nframes = None

    # data_param = f"data/{phi}_EKRM_trajs.dat"
    # window_size_x=217.6
    # window_size_y=174
    # Nframes=2400

    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles  = data['particles']
    time_step  = data['time_step']
    pixel_size = data['pixel_size']
    particles[:, [2]] += 1 # make t 1-based
    Nframes = None
    window_size_x = None
    window_size_y = None

    
    box_sizes = np.array([1, 2, 4, 8, 16, 32, 64, 128]) * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    sep_sizes = np.full_like(box_sizes, 2*a)
    # sep[0] = -60

    N2_mean, N2_std, N_stats = countoscope.Calc_and_Output_Stats(data=particles, 
                                                                 Nframes=Nframes, 
                                                                 window_size_x=window_size_x, window_size_y=window_size_y, 
                                                                 box_sizes=box_sizes,
                                                                #  box_sizes_x=box_sizes_x, box_sizes_y=box_sizes_y,
                                                                 sep_sizes=sep_sizes,)
    np.savez(f'box_counting/data/counted_{file}.npz', N2_mean=N2_mean, N2_std=N2_std, N_stats=N_stats,
             box_sizes=box_sizes, sep_sizes=sep_sizes,
             time_step=time_step)
    
t1 = time.time()
print(f'done in {t1-t0:.0f}s')