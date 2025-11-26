import numpy as np
import common

if __name__ == '__main__':
    data = np.loadtxt(f'raw_data/tristan/V_sweep_9V_10khz_014_tracks.txt')
    pixel_size = 0.11

    data[:, 2] -= data[:, 2].min() # ensure time is zero based

    data[:, [0, 1]] *= pixel_size

    common.save_data(f'particle_detection/data/particles_tristan0.npz', particles=data,
                time_step=1, num_timesteps=data[:, 2].max()+1,
                pixel_size=pixel_size,
                window_size_x=data[:, 0].max(), window_size_y=data[:, 1].max())