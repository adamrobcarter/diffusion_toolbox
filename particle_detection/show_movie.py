import numpy as np
import common
import sys
import particle_detection.show
import preprocessing.stack_movie

if __name__ == '__main__':
    # SUFFIX = '_nominmass'
    SUFFIX = ''

    for file in sys.argv[1:]:

        data_particles = common.load(f'particle_detection/data/particles_{file}{SUFFIX}.npz')
        particles = data_particles['particles']
        radius    = data_particles.get('radius')
        time_step = data_particles['time_step']
        num_timesteps = particles[:, 2].max() + 1
        window_size_x = data_particles['window_size_x']
        window_size_y = data_particles['window_size_y']

        try:
            data_stack = common.load(f'preprocessing/data/stack_{file}.npz')

        except FileNotFoundError:
            num_timesteps = int(particles[:, 2].max() - 1)
            stack = None
            pixel_size = 1
            radius = np.full(particles.shape[0], 0.002*160)

            no_stack = True
        
        else:
            no_stack = False

            stack = data_stack['stack']
            pixel_size = data_stack['pixel_size']

            # crop
            stack = stack[:, :500, :500]

            # stack = common.add_drift_intensity(stack, 1)

            print(stack.shape[1], 'x', stack.shape[2], 'px')


            stack = stack - stack.mean(axis=0) # remove space background
        
        def add_outlines(timestep, ax):
            particle_detection.show.add_particle_outlines(ax, pixel_size, particles, radius, timestep)

        filename = f'movie_{file}{SUFFIX}.gif'
        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"preprocessing/figures_png/{filename}",
                                                func=add_outlines, num_timesteps_in_data=num_timesteps,
                                                window_size_x=window_size_x, window_size_y=window_size_y)
        # preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}{SUFFIX}.mp4",
        #                                            func=add_outlines)
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/{filename}.mp4")