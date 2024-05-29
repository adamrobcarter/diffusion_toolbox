import numpy as np

data = np.loadtxt(f'raw_data/tristan/V_sweep_9V_10khz_014_tracks.txt')

data[:, 2] -= data[:, 2].min() # ensure time is zero based

np.savez(f'particle_detection/data/particles_tristan0.npz', particles=data,
            time_step=1, num_timesteps=data[:, 2].max()+1)