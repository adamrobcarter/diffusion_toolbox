import numpy as np
import scattering_functions.scattering_functions_nonumba as scattering_functions
import common
import time

log = True

# num_k_bins = 200 # was 50
# num_k_bins = int(max_K / min_K)
# num_k_bins = 100
num_k_bins = 100
num_iters = 24 # was 10
#                                       9.1 so that 9 gets included                round as we need integer frames

drift_removed = False
# crop = 0.5 if F_type == 'F' else 1.0
# crop = 0.5 # to force the same
crop = 1.0

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
       
        num_timesteps = particles[:, 2].max()
        # d_frames = np.concatenate([np.arange(0, 9.1), np.logspace(1, np.log10(num_timesteps-100), 50)]).round()
        d_frames = common.exponential_integers(1, num_timesteps-1) - 1 

        if drift_removed:
            particles = common.remove_drift(particles)


        if False:
            print('adding drift')
            particles = common.add_drift(particles, 0.05, 0.05)

        width  = particles[:, 0].max() - particles[:, 0].min() # what are these for?
        height = particles[:, 1].max() - particles[:, 1].min()


        # min_K = 2 * np.pi / (min(width, height) * 0.5) # 0.5 was crop but then min_K varies with F_type
        # min_K = min_K * 2
        # print("min_k", min_K)
        #print(f"min K = {min_K:.3f}")
        max_K = 10 # was 10

        # max_time_origins = 50
        max_time_origins = 100
        max_time_origins = 1000

        print('starting calculation')

        Fs, F_unc, ks = scattering_functions.intermediate_scattering(log, F_type, crop, num_k_bins, max_time_origins, d_frames, 
                                                                     particles, max_K, width=width, height=height)

        t1 = time.time()
            
        common.save_data(f"scattering_functions/data/{F_type}_{file}", F=Fs, F_unc=F_unc, k=ks, t=d_frames*time_step, crop=crop,
                num_k_bins=num_k_bins, max_time_origins=max_time_origins, computation_time=t1-t0, log=log,
                particle_diameter=particle_diameter, drift_removed=drift_removed)

        print()

if __name__ == '__main__':
    calc_for_f_type('F_s')
    calc_for_f_type('f')
