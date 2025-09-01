import numpy as np
import tqdm

def calc(particles, grid_size, window_size_x, window_size_y):
    # sort by ID then time
    particles = particles[np.lexsort((particles[:, 2], particles[:, 3]))]
    
    # for row_i in tqdm.trange(particles.shape[0]):
    num_timesteps = int(particles[:, 2].max() + 1)

    progress = tqdm.tqdm(total=particles.shape[0])

    grid_xs = np.arange(0, window_size_x, grid_size)
    grid_ys = np.arange(0, window_size_y, grid_size)

    v_x_sum = np.zeros((grid_xs.size, grid_ys.size, num_timesteps))
    v_y_sum = np.zeros((grid_xs.size, grid_ys.size, num_timesteps))
    n       = np.zeros((grid_xs.size, grid_ys.size, num_timesteps))

    row_i = 0
    while row_i < particles.shape[0]:
        this_id = particles[row_i, 3]

        while row_i < particles.shape[0]-1 and particles[row_i+1, 3] == this_id:
            dx = particles[row_i+1, 0] - particles[row_i, 0]
            dy = particles[row_i+1, 1] - particles[row_i, 1]
            centre_x = (particles[row_i+1, 0] + particles[row_i, 0])/2
            centre_y = (particles[row_i+1, 1] + particles[row_i, 1])/2
            index_x = np.argmax(grid_xs > centre_x)
            index_y = np.argmax(grid_ys > centre_y)
            v_x_sum[index_x, index_y, int(particles[row_i, 2])] += dx
            v_y_sum[index_x, index_y, int(particles[row_i, 2])] += dy
            n      [index_x, index_y, int(particles[row_i, 2])] += 1

            row_i += 1
            progress.update()
        
        # this shouldn't be here but something needs to I think in case the while evaluates false on the first try
        row_i += 1
        progress.update()


    return v_x_sum/n, v_y_sum/n, grid_xs, grid_ys, n

