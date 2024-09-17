import numpy as np
import scattering_functions.scattering_functions_nonumba as scattering_functions
import common
import time
import warnings

log = True
if not log:
    warnings.warn('not using log calculation')

# num_k_bins = 100
num_k_bins = 50
# computation time is proportional to this squared

drift_removed = False



max_time_origins = 1000
# max_time_origins = 100
# this sets (approx) how many different time origins will be averaged over
# computation time is directly proportional

def calc_for_f_type(F_type):

    for file in common.files_from_argv('particle_detection/data', 'particles_'):

        t0 = time.time()

        if F_type == 'Fs':
            filepath = f"particle_linking/data/trajs_{file}.npz"
        else:
            filepath = f"particle_detection/data/particles_{file}.npz"
        data = common.load(filepath)
        particles         = data['particles'] # rows of x,y,t
        time_step         = data['time_step']
        particle_diameter = data.get('particle_diameter')
        pixel_size=data.get('pixel_size')
       
        # num_timesteps = particles[:, 2].max()
        num_timesteps = data['num_timesteps']

        # d_frames = np.concatenate([np.arange(0, 9.1), np.logspace(1, np.log10(num_timesteps-100), 50)]).round()
        d_frames = common.exponential_integers(1, num_timesteps-1) - 1 

        if drift_removed:
            particles = common.remove_drift(particles)
        else:
            print('not removing drift')

        if False:
            print('adding drift')
            particles = common.add_drift(particles, 0.05, 0.05)

        width  = particles[:, 0].max() - particles[:, 0].min() # what are these for?
        height = particles[:, 1].max() - particles[:, 1].min()
        
        min_K = 2*np.pi/min(width, height)
        if file.startswith('brennan') or file == 'eleanorlong':
            min_K = 2*np.pi/293.216291078
            # we use this cause it's the same as used by eleanorlong

        if pixel_size:
            max_K = 2 * np.pi / pixel_size
        else:
            print('Pixel size not given, using max_K = 21.81661564992912')
            max_K = 21.81661564992912 # was 10
            # we use this cause it's the same as used by eleanorlong

        print('starting calculation')

        Fs, F_unc, ks = scattering_functions.intermediate_scattering(log, F_type, num_k_bins, max_time_origins, d_frames, 
                                                                     particles, max_K, min_K)

        t1 = time.time()
            
        filename = f"scattering_functions/data/{F_type}_{file}"
        if not log:
            filename += '_nolog'
        common.save_data(filename, F=Fs, F_unc=F_unc, k=ks, t=d_frames*time_step,
            num_k_bins=num_k_bins, max_time_origins=max_time_origins, computation_time=t1-t0, log=log,
            particle_diameter=particle_diameter, drift_removed=drift_removed,
            pixel_size=pixel_size, pack_frac_given=data.get('pack_frac_given'),
            window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
            NAME=data.get('NAME'), channel=data.get('channel'),
        )

        print()

if __name__ == '__main__':
    calc_for_f_type('F_s')
    calc_for_f_type('f')
