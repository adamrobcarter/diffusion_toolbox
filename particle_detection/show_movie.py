import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.animation
import tqdm
import particle_detection.show


# def go(filename, stack_file, particles_file, output_file):
#     # this is also used by particle_linking.show_movie
    
#     data = common.load(stack_file)
#     stack             = data['stack']
#     pixel_size        = data['pixel_size']
#     # particle_diameter = data['particle_diameter']
#     # num_timesteps = stack.shape[0]

#     # print(stack.shape[1], 'x', stack.shape[2], 'px')

#     stack = stack - stack.mean(axis=0)
#     stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

#     # print('stack min max', stack.min(), stack.max())

#     data = common.load(particles_file)
#     particles = data['particles']
#     radius    = data['radius']
    
#     print(f'found {(particles.shape[0]/particles[:, 2].max()):0f} particles per frame')
#     # print(particles[:, 2].min(), particles[:, 2].max())

#     fig, ax = plt.subplots(1, 1)

#     end_frame = min(stack.shape[0], 50)
#     frames = range(0, end_frame, 1)

#     def show(timestep):
#         ax.clear()

#         particle_detection.show.show_frame(fig, ax, stack, pixel_size, particles, radius, timestep, filename)
#         ax.set_title(f'{filename} t={timestep:03d}')

#     common.save_gif(show, frames, fig, output_file, fps=2)

# if __name__ == '__main__':
#     for file in common.files_from_argv('preprocessing/data', 'stack_'):
#         go(file,
#            stack_file=f'preprocessing/data/stack_{file}.npz',
#            particles_file=f'particle_detection/data/particles_{file}.npz',
#            output_file=f"particle_detection/figures_png/movie_{file}.gif")

import preprocessing.stack_movie
import particle_detection.show


if __name__ == '__main__':
    # SUFFIX = '_nominmass'
    SUFFIX = ''

    for file in sys.argv[1:]:
        data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
        
        stack = data_stack['stack']
        pixel_size = data_stack['pixel_size']
        time_step = data_stack['time_step']

        
        data_particles = common.load(f'particle_detection/data/particles_{file}{SUFFIX}.npz')
        particles = data_particles['particles']
        radius    = data_particles['radius']

        # crop
        stack = stack[:, :500, :500]

        # stack = common.add_drift_intensity(stack, 1)

        print(stack.shape[1], 'x', stack.shape[2], 'px')


        stack = stack - stack.mean(axis=0) # remove space background
        
        def add_outlines(timestep, ax):
            particle_detection.show.add_particle_outlines(ax, pixel_size, particles, radius, timestep)

        filename = f'movie_{file}'
        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"preprocessing/figures_png/{filename}{SUFFIX}.gif",
                                                func=add_outlines)
        # preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}{SUFFIX}.mp4",
        #                                            func=add_outlines)
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/{filename}.mp4")