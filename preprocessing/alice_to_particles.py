import numpy as np

for phi in [0.02, 0.34, 0.66]:
    data = np.fromfile(f'raw_data/{phi}_newtrack.dat', dtype=float, sep=' ') # Alice
    #all_data = np.loadtxt('data/0.34_EKRM_trajs.dat', delimiter=',', skiprows=1) # Eleanor
    data = data.reshape((-1, 4))
    PIXEL = 0.17
    data[:, [0,1]] *= PIXEL

    num_timesteps = data[:, 2].max() - data[:, 2].min()

    data[:, 2] -= data[:, 2].min() # ensure time is zero based
    data[:, 3] -= data[:, 3].min() # ensure ID is zero based

    np.savez(f'particle_detection/data/particles_alice{phi}.npz', particles=data,
             time_step=0.5, particle_diameter=2.8, pixel_size=PIXEL,
             num_timesteps=num_timesteps)