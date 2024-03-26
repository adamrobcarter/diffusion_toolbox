import countoscope
import numpy as np
import time
import sys
import common

def calc_and_save(box_sizes_px, sep_sizes_px, file, output_file_name):
    t0 = time.time()

    a = 1.395 #radius of particles

    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles                  = data['particles']
    time_step                  = data['time_step']
    pixel_size                 = data['pixel_size']
    num_timesteps              = data['num_timesteps']
    particle_diameter          = data['particle_diameter']
    # particle_diameter_calced   = data['particle_diameter_calced']
    particle_diameter_calced   = data.get('particle_diameter_calced')
    # particles[:, [2]] += 1 # make t 1-based
    Nframes = None
    window_size_x = None
    window_size_y = None

    if file == 'eleanorlong':
        pixel_size = 0.17 # so the boxes are the same size as aliceXXX
    
    box_sizes = box_sizes_px * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    sep_sizes = sep_sizes_px * pixel_size 
    # bigger sep at low box sizes to reduce the number of boxes, as we already have loads of stats for small boxes
    # sep[0] = -60
    
    added_drift_x = 0
    added_drift_y = 0
    if False:
        added_drift_x = 0.03
        particles = common.add_drift(particles, added_drift_x, added_drift_y)

    N2_mean, N2_std, N_stats = countoscope.calculate_N1N2(data=particles,
                                                          box_sizes=box_sizes)
    

    # also let's calc the packing fraction now as it'll be useful in the future
    window_width  = particles[:, 0].max() - particles[:, 0].min()
    window_height = particles[:, 1].max() - particles[:, 1].min()
    density = particles.shape[0]/num_timesteps / (window_width * window_height)
    pack_frac = np.pi/4 * density * particle_diameter**2

    np.savez(output_file_name, N2_mean=N2_mean, N2_std=N2_std,
             N_stats=N_stats, box_sizes=box_sizes, sep_sizes=sep_sizes,
             time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
             added_drift_x=added_drift_x, added_drift_y=added_drift_y,
             particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0)
    # np.savez(f'box_counting/data/all_counts_{file}.npz', counts=counts,
    #          N_stats=N_stats, box_sizes=box_sizes, sep_sizes=sep_sizes,
    #          time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
    #          particle_diameter_calced=particle_diameter_calced)

if __name__ == '__main__':
    for file in sys.argv[1:]:
        # box_sizes_px = np.array([ 1,  2,  4,  8, 16, 32, 64, 128])
        # # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
        # sep_sizes_px = np.array([50, 50, 50, 50, 50, 50, 50,  50])
        filename = f'box_counting/data/countedN1N2_{file}_nodrift.npz'

        box_sizes_px = np.array([ 2,  4,  8, 16, 32])
        sep_sizes_px = np.array([20, 20, 20, 20, 20])
        calc_and_save(box_sizes_px, sep_sizes_px, file, filename)

        # box_sizes_px = np.array([16, 32, 64, 128])
        # # sep_sizes_px = np.array([ 0,  0,  0,  0])
        # sep_sizes_px = -box_sizes_px/2
        # filename = f'box_counting/data/counted_{file}_overlapped_neg.npz'
        # calc_and_save(box_sizes_px, sep_sizes_px, file, filename)

        # box_sizes_px = np.array([16, 32, 64, 128])
        # sep_sizes_px = np.array([-10, -20, -20, -20])
        # filename = f'box_counting/data/counted_{file}_overlapped2.npz'
        # calc_and_save(box_sizes_px, sep_sizes_px, file, filename)

        # box_sizes_px = np.array([16, 32, 64, 128])
        # sep_sizes_px = np.array([-10, -25, -40, -100])
        # filename = f'box_counting/data/counted_{file}_overlapped3.npz'
        # calc_and_save(box_sizes_px, sep_sizes_px, file, filename)