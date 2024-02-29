import numpy as np

for phi in [0.02, 0.34, 0.66]:
    data_param = np.fromfile(f'raw_data/{phi}_newtrack.dat', dtype=float, sep=' ') # Alice
    #all_data = np.loadtxt('data/0.34_EKRM_trajs.dat', delimiter=',', skiprows=1) # Eleanor
    data_param = data_param.reshape((-1, 4))
    PIXEL = 0.17
    data_param[:, [0,1]] *= PIXEL

    num_timesteps = data_param[:, 2].max() - data_param[:, 2].min()

    np.savez(f'particle_detection/data/particles_alice{phi}.npz', particles=data_param,
             time_step=0.5, particle_diameter=2.8, pixel_size=PIXEL,
             num_timesteps=num_timesteps)