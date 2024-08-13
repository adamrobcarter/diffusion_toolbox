import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.animation
import tqdm
import particle_detection.show
import matplotlib.cm

import preprocessing.stack_movie
import particle_detection.show


if __name__ == '__main__':
    # SUFFIX = '_nominmass'
    SUFFIX = ''

    for file in sys.argv[1:]:

        # data_particles = common.load(f'particle_detection/data/particles_{file}{SUFFIX}.npz')
        # particles = data_particles['particles']
        # radius    = data_particles.get('radius')
        # time_step = data_particles['time_step']

        data_counts = common.load(f'box_counting/data/counted_{file}_for_equilib.npz')
        L = data_counts['box_sizes'][0]
        avg_density = data_counts['N_mean'][0] / L**2
        box_coords = data_counts['box_coords'][0, :, :, :]
        xs = box_coords[:, :, 0]
        ys = box_coords[:, :, 1]
        time_step = data_counts['time_step']

        try:
            data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
        
            stack = data_stack['stack']
            pixel_size = data_stack['pixel_size']

            # crop
            stack = stack[:, :500, :500]

            # stack = common.add_drift_intensity(stack, 1)

            print(stack.shape[1], 'x', stack.shape[2], 'px')


            stack = stack - stack.mean(axis=0) # remove space background

        except FileNotFoundError:
            pass
            # num_timesteps = int(particles[:, 2].max() - 1)
            # stack = np.zeros((num_timesteps, 320, 320))
            # pixel_size = 1
            # radius = np.full(particles.shape[0], 0.002*160)
        
        def add_outlines(timestep, ax):
            # particle_detection.show.add_particle_outlines(ax, pixel_size, particles, radius, timestep)
            counts = data_counts['counts'][0, :, :, timestep]
            for x in range(counts.shape[0]):
                for y in range(counts.shape[1]):
                    density = counts[x, y] / L**2
                    relative_density = density / avg_density
                    color = matplotlib.cm.RdBu(np.interp(relative_density, (0.95, 1.05), (0, 1)))
                    rect = matplotlib.patches.Rectangle((xs[y, x], ys[y, x]), L, L,
                        linewidth=1, facecolor=color, alpha=0.4)
                    ax.add_patch(rect)

        filename = f'movie_{file}'
        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"box_counting/figures_png/equilib_{filename}.gif",
                                                func=add_outlines)
        # preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}{SUFFIX}.mp4",
        #                                            func=add_outlines)
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/{filename}.mp4")