import numpy as np
import scattering_functions.scattering_functions_nonumba as scattering_functions
import common
import time
import warnings

def calc_for_f_type(
        file,
        F_type,
        log=True,
        max_K=None,
        drift_removed=False,
        num_k_bins=50, # computation time is proportional to this squared
        file_suffix='',
        cores=16,
        max_time_origins=100, # this sets (approx) how many different time origins will be averaged over
                              # computation time is directly proportional
        S_only=False,
        first_point_only=False,
        use_zero=False,
        save=True,
    ):


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
        width             = data['window_size_x']
        height            = data['window_size_y']
       
        num_timesteps = int(particles[:, 2].max() + 1)

        # d_frames = np.concatenate([np.arange(0, 9.1), np.logspace(1, np.log10(num_timesteps-100), 50)]).round()
        if S_only:
            d_frames = np.array([0])
        elif first_point_only:
            d_frames = np.array([0, 1])
        else:
            d_frames = common.exponential_integers(1, num_timesteps-1) - 1 

        if drift_removed:
            particles = common.remove_drift(particles)
        else:
            print('not removing drift')

        if False:
            print('adding drift')
            particles = common.add_drift(particles, 0.05, 0.05)
        
        min_K = 2*np.pi/min(width, height)
        if file.startswith('brennan'):
            min_K = 2*np.pi/293.216291078
            # we use this cause it's the same as used by eleanorlong
        print(f'min k = {min_K:.3f}')

        if not max_K:
            if pixel_size:
                max_K = 2 * np.pi / pixel_size
            else:
                print('Pixel size not given, using max_K = 21.81661564992912')
                max_K = 21.81661564992912 # was 10
                # we use this cause it's the same as used by eleanorlong

        print('starting calculation')

        Fs, F_unc, ks, F_unbinned, F_unc_unbinned, k_unbinned, k_x, k_y = scattering_functions.intermediate_scattering(log, F_type, num_k_bins, max_time_origins, d_frames, 
                                                                     particles, max_K, min_K, cores=cores, use_zero=use_zero)

        t1 = time.time()
            
        filename = f"scattering_functions/data/{F_type}_{file}{file_suffix}"
        # if not log:
        #     filename += '_nolog'

        to_save = dict(F=Fs, F_unc=F_unc, k=ks, k_x=k_x, k_y=k_y,
                F_unbinned=F_unbinned, k_unbinned=k_unbinned, F_unc_unbinned=F_unc_unbinned,
                t=d_frames*time_step,
                num_k_bins=num_k_bins, max_time_origins=max_time_origins, computation_time=t1-t0, log=log,
                particle_diameter=particle_diameter, drift_removed=drift_removed,
                pixel_size=pixel_size, pack_frac_given=data.get('pack_frac_given'),
                window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
                NAME=data.get('NAME'), channel=data.get('channel'),
            )

        if save:
            common.save_data(filename, **to_save)

        print()

        return to_save

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        calc_for_f_type(file, 'F_s')
        calc_for_f_type(file, 'f')
