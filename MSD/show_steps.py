import common
import numpy as np
import matplotlib.pyplot as plt
import tqdm

for file in common.files_from_argv('particle_linking/data', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']
    
    particles = particles[np.lexsort((particles[:, 2], particles[:, 3]))] # sort by ID then time

    x_steps = np.full(particles.shape[0], np.nan)
    y_steps = np.full(particles.shape[0], np.nan)
    num_steps = 0

    row_i = 0
    progress = tqdm.tqdm(total=particles.shape[0])

    while row_i < particles.shape[0]:
        current_id = particles[row_i, 3]

        while row_i < particles.shape[0]-1 and particles[row_i+1, 3] == current_id:
            row_i += 1
            progress.update()

            dr = particles[row_i, [0, 1]] - particles[row_i-1, [0, 1]]
            x_steps[num_steps] = dr[0]
            y_steps[num_steps] = dr[1]
            num_steps += 1
        else:
            row_i += 1
            progress.update()
        
    progress.close()
    
    x_steps = x_steps[np.isfinite(x_steps)]
    y_steps = y_steps[np.isfinite(y_steps)]

    fig, (ax_x, ax_y) = plt.subplots(2, 1, figsize=(4, 6))
    max_step = max(np.abs(x_steps).max(), np.abs(y_steps).max())
    bins = np.linspace(-max_step, max_step, 21)
    ax_x.hist(x_steps, bins=bins)
    ax_y.hist(y_steps, bins=bins)
    ax_x.semilogy()
    ax_y.semilogy()
    common.save_fig(fig, f'MSD/figures_png/steps_{file}.png')

