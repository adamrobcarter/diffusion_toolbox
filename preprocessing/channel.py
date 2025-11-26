import numpy as np
import common

if __name__ == '__main__':
    # source_file = 'raw_data/eleanor/20fps_s2.txt'
    source_file = 'raw_data/eleanor/30fps_s6.txt'
    data = np.loadtxt(source_file, skiprows=1)

    # xy are flipped
    # data = data[:, [1, 0, 2]]

    # time starts at zero
    data[:, 2] -= data[:, 2].min()
    last_timestep = data[:, 2].max()

    pixel_size = 0.3
    data[:, [0, 1]] *= pixel_size
    dt = 1/30
    window_size_x = 1280 * pixel_size
    window_size_y = 1024 * pixel_size
    particle_diameter = 2.8
    outfile = 'channel'

    assert np.all(data[:, 0] < window_size_x)
    assert np.all(data[:, 1] < window_size_y)

    print(data.shape)
    for i in [0, 1, 2]:
        print(data[:, i].min(), data[:, i].max())

    common.save_data(f'particle_detection/data/particles_{outfile}.npz',
        particles=data,
        pixel_size=pixel_size,
        time_step=dt, particle_diameter=particle_diameter,
        window_size_x=window_size_x, window_size_y=window_size_y, max_time_hours=round(last_timestep*dt/60/60, 2),
        source_file=source_file, 
        #density=density, extra_source_file=extra_source_file,
    )