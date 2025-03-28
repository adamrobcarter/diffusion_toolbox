import numpy as np
import common
import sys
import particle_detection.show
import preprocessing.stack_movie

CROP = None
HIGHLIGHTS = True

def go(file, infile, outfile, SUFFIX=''):

        data_particles = common.load(infile)
        particles = data_particles['particles']
        # particles = particles[:, [1, 0, 2]] # DON'T ASK ME WHY

        print('particles max', particles[:, 0].max(), particles[:, 0].min())
        print('particles max', particles[:, 1].max(), particles[:, 1].min())
        
        radius    = data_particles.get('radius')
        time_step = data_particles['time_step']
        num_timesteps = particles[:, 2].max() + 1
        window_size_x = data_particles['window_size_x']
        window_size_y = data_particles['window_size_y']

        if CROP:
            in_crop = (particles[:, 0] < CROP) & (particles[:, 1] < CROP)
            particles = particles[in_crop, :]
            window_size_x = CROP
            window_size_y = CROP

        try:
            data_stack = common.load(f'preprocessing/data/stack_{file}.npz')

        except FileNotFoundError:
            num_timesteps = int(particles[:, 2].max() - 1)
            stack = None
            pixel_size = None
            # radius = np.full(particles.shape[0], 0.002*160)

            # no_stack = True
        
        else:
            # no_stack = False

            stack = data_stack['stack']
            pixel_size = data_stack['pixel_size']

            
            stack = stack[:, ::-1, :]

            # crop
            if CROP:
                stack = stack[:, :CROP*pixel_size, :500*pixel_size]

            # stack = common.add_drift_intensity(stack, 1)

            print(stack.shape[1], 'x', stack.shape[2], 'px')


            stack = stack - stack.mean(axis=0) # remove space background
            
        if 'particle_diameter' in data_particles:
            radius = np.full(particles.shape[0], data_particles['particle_diameter']/2)
        else:
            radius = None
        
        def add_outlines(timestep, ax):
            particle_detection.show.add_particle_outlines(ax, pixel_size, particles, radius, timestep, outline=False)

        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, outfile,
                                                func=add_outlines, num_timesteps_in_data=num_timesteps,
                                                window_size_x=window_size_x, window_size_y=window_size_y,
                                                highlights=HIGHLIGHTS)
        # preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}{SUFFIX}.mp4",
        #                                            func=add_outlines)
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/{filename}.mp4")



if __name__ == '__main__':

    for file in sys.argv[1:]:
        go(
            file,
            infile = f'particle_detection/data/particles_{file}.npz',
            outfile = f'particle_detection/figures_png/movie_{file}.gif'
        )