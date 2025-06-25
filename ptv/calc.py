import common
import ptv.ptv
import numpy as np

for file in common.files_from_argv('particle_linking/data', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']
    time_step = data['time_step']
    window_size_x = data['window_size_x']
    window_size_y = data['window_size_y']

    grid_size = 40
    v, grid_xs, grid_ys = ptv.ptv.calc(particles, grid_size=grid_size, window_size_x=window_size_x, window_size_y=window_size_y)

    v = np.nanmean(v, axis=2) # average over time

    common.save_data(f'ptv/data/ptv_{file}.npz', v=v, time_step=time_step, grid_size=grid_size,
                     grid_xs=grid_xs, grid_ys=grid_ys)