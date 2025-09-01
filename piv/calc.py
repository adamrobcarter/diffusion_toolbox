import common
import piv.piv
import numpy as np

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    time_step = data['time_step']

    particle_size_px = int(data['particle_diameter'] / data['pixel_size'])
    integration_window = particle_size_px*2
    print(f'interogation_window = {integration_window}px')
    v_x, v_y, signal_to_noise = piv.piv.go(stack, time_step, interrogation_window=integration_window)

    print(v_x.shape, v_y.shape)

    # v_x = np.nanmean(v_x, axis=0)
    # v_y = np.nanmean(v_y, axis=0)
    n = np.isfinite(v_x).sum(axis=0) # assume v_x and v_y are nan in the same places?

    common.save_data(f'piv/data/piv_{file}.npz',
        v_x=v_x, v_y=v_y, signal_to_noise=signal_to_noise,
        n=n, time_step=time_step,
        integration_window=integration_window,
    #  grid_size=grid_size,
        # grid_xs=grid_xs, grid_ys=grid_ys,
        particle_diameter=data.get('particle_diameter'), particle_material=data.get('particle_material'),
        # window_size_x=data['window_size_x'], window_size_y=data['window_size_y']
    )