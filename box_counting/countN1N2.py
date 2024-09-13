import countoscope_old as countoscope
import numpy as np
import time
import sys
import common

def calc_and_save(box_sizes_px, sep_sizes_px, file, output_file_name, drift):
    t0 = time.time()

    a = 1.395 #radius of particles

    # data = common.load(f'particle_detection/data/particles_{file}.npz')
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles                  = data['particles']
    time_step                  = data['time_step']
    pixel_size                 = data['pixel_size']
    num_timesteps              = data['num_timesteps']
    particle_diameter          = data['particle_diameter']
    particle_diameter_calced   = data.get('particle_diameter_calced')
    window_size_x            = data.get('window_size_x')
    window_size_y            = data.get('window_size_y')

    if file == 'eleanorlong':
        pixel_size = 0.17 # so the boxes are the same size as aliceXXX
    
    box_sizes = box_sizes_px * pixel_size 
    sep_sizes = sep_sizes_px * pixel_size 
    # ^^ there's no reason for the boxes to be integer pixel multiples, but it helps comparisons with the intensity method
    
    t0 = time.time()
    particles = common.remove_drift(particles) # first make sure there's no drift
    t1 = time.time()
    print('took', t1-t0)

    added_drift_x = 0
    added_drift_y = 0
    if drift != 0:
        added_drift_x = drift
        particles = common.add_drift(particles, added_drift_x, added_drift_y)
        
    drift_x = added_drift_x / time_step # go from um/frame to um/s
    drift_y = added_drift_y / time_step

    results = countoscope.calculate_N1N2(data=particles, box_sizes=box_sizes, sep_sizes=sep_sizes,
                                        #  window_size_x=window_size_x, window_size_y=window_size_y
                                         )
    
    # also let's calc the packing fraction now as it'll be useful in the future
    window_width  = particles[:, 0].max() - particles[:, 0].min()
    window_height = particles[:, 1].max() - particles[:, 1].min()
    density = particles.shape[0]/num_timesteps / (window_width * window_height)
    pack_frac = np.pi/4 * density * particle_diameter**2

    common.save_data(output_file_name, N1N2_mean=results.N1N2, N1N2_std=results.N1N2_std,
             N_mean=results.N_mean, N_var=results.N_var,
             box_coords_1=results.box_coords_1, box_coords_2=results.box_coords_2,
             box_sizes=box_sizes, sep_sizes=sep_sizes,
             time_step=time_step, pack_frac=pack_frac, particle_diameter=particle_diameter,
             drift_x=drift_x, drift_y=drift_y,
             particle_diameter_calced=particle_diameter_calced, computation_time=time.time()-t0,
             window_size_x=window_size_x, window_size_y=window_size_y)


if __name__ == '__main__':
    for file in sys.argv[1:]:
        # box_sizes_px = np.array([32,  64,  96,  128])
        box_sizes_sigma = np.array([0.25, 0.5, 1, 2, 4, 8])
        box_sizes_px = box_sizes_sigma * 2.8 / 0.17
        sep_sizes_px = 40 - box_sizes_px
        
        # box_sizes_px = np.array([4,   8, 16, 22, 32])
        # sep_sizes_px = np.array([20, 15, 5, 5, 5])

        filename = f'box_counting/data/countedN1N2_{file}_nodrift.npz'
        calc_and_save(box_sizes_px, sep_sizes_px, file, filename, 0.0)

        filename = f'box_counting/data/countedN1N2_{file}_drift0.03.npz'
        calc_and_save(box_sizes_px, sep_sizes_px, file, filename, 0.05)