import common
import ptv.ptv
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particles = data['particles']
        time_step = data['time_step']
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']

        grid_size = 10
        v_x, v_y, grid_xs, grid_ys, n = ptv.ptv.calc(particles, grid_size=grid_size, window_size_x=window_size_x, window_size_y=window_size_y)

        v_x = np.nanmean(v_x, axis=2) # average over time
        v_y = np.nanmean(v_y, axis=2) # average over time
        n   = np.nanmean(n,   axis=2) # average over time
        v_x *= data['pixel_size'] / data['time_step']
        v_y *= data['pixel_size'] / data['time_step']

        common.save_data(f'ptv/data/ptv_{file}.npz', v_x=v_x, v_y=v_y, n=n, time_step=time_step, grid_size=grid_size,
            grid_xs=grid_xs, grid_ys=grid_ys,
            particle_diameter=data.get('particle_diameter'), particle_material=data.get('particle_material'),
            window_size_x=data['window_size_x'], window_size_y=data['window_size_y']
        )