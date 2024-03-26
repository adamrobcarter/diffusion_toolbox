import numpy as np
import scattering_functions.scattering_functions_nonumba as scattering_functions
import common
import time
import sys

log = True

# num_k_bins = 200 # was 50
# num_k_bins = int(max_K / min_K)
# num_k_bins = 100
num_k_bins = 100
num_iters = 24 # was 10
#                                       9.1 so that 9 gets included                round as we need integer frames

F_types = ['F', 'F_s']
F_types = ['F']

drift_removed = False

for file in sys.argv[1:]:
    # crop = 0.5 if F_type == 'F' else 1.0
    # crop = 0.5 # to force the same
    crop = 1.0

    t0 = time.time()

    data = common.load(f"particle_detection/data/particles_{file}.npz")
    particles = data['particles']
    time_step = data['time_step']
    particle_diameter = data['particle_diameter']
    # rows of x,y,t
    print(particles.shape)
    num_timesteps = particles[:, 2].max()
    d_frames = np.concatenate([np.arange(0, 9.1), np.logspace(1, np.log10(num_timesteps-100), 50)]).round()
    # d_frames = np.concatenate([np.arange(0, 9.1), np.logspace(1, np.log10(1500), 50)]).round()


    if drift_removed:
        particles = common.remove_drift(particles)

    # if drift_removed:
    #     data = remove_drift(data)

    # num_timesteps = particles.shape[1]

    # width  = data['width']
    # height = data['height']
    # num_frames = particles.shape[1]

    particles = common.add_drift(particles, 0.05, 0.05)

    width  = particles[:, 0].max() - particles[:, 0].min() # what are these for?
    height = particles[:, 1].max() - particles[:, 1].min()

    for F_type in F_types:

    #     print(particles.shape)
    #     particles[:, 0] -= particles[:, 0].mean()
    #     particles[:, 1] -= particles[:, 1].mean()
    #     print(particles[:, 0].mean(), 'vs', width/2)
    #     print(particles[:, 1].mean(), 'vs', height/2)
    #     print(particles[:, 0].mean())
    #     print(particles[:, 1].mean())
    #     print(particles.shape)

        # min_K = 2 * np.pi / (min(width, height) * 0.5) # 0.5 was crop but then min_K varies with F_type
        # min_K = min_K * 2
        # print("min_k", min_K)
        #print(f"min K = {min_K:.3f}")
        max_K = 10 # was 10

        # values of d_frame to use. Don't need as much resolution at large t, but need integers,
        # hence this slightly weird thing. Could the whole thing just be a logspace?
        #d_frames = [0]

        Fs, F_unc, ks = scattering_functions.intermediate_scattering(log, F_type, crop, num_k_bins, num_iters, d_frames, particles, max_K, width=width, height=height)

        # suffix = '_loglog' if log else ''

        t1 = time.time()
            
        np.savez(f"scattering_functions/data/{F_type}_{file}", F=Fs, F_unc=F_unc, k=ks, t=d_frames*time_step, crop=crop,
                num_k_bins=num_k_bins, num_iters=num_iters, computation_time=t1-t0, log=log,
                particle_diameter=particle_diameter, drift_removed=drift_removed)

        print()
