import numpy as np
import common
import sys
import particle_detection.show
import preprocessing.stack_movie

CROP = None
HIGHLIGHTS = False

def go(file, infile, outfile, SUFFIX='', crop=False, every_nth_frame=None, output_type='movie',
       annotation_color='white',):

        data_particles = common.load(infile)
        particles = data_particles['particles']
        # particles = particles[:, [1, 0, 2]] # DON'T ASK ME WHY

        print('particles max', particles[:, 0].max(), particles[:, 0].min())
        print('particles max', particles[:, 1].max(), particles[:, 1].min())
        
        time_step = data_particles['time_step']
        num_timesteps = particles[:, 2].max() + 1
        assert num_timesteps > 0, f'num_timesteps = {num_timesteps}'
        window_size_x = data_particles['window_size_x']
        window_size_y = data_particles['window_size_y']

        try:
            data_stack = common.load(f'preprocessing/data/stack_{file}.npz')

        except FileNotFoundError:
            num_timesteps = int(particles[:, data_particles.get('dimension', 2)].max() - 1)
            assert num_timesteps > 0, f'num_timesteps = {num_timesteps}'
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
            if crop:
                print('please crop the stack properly')
                # stack = stack[:, int((window_size_x - CROP)/pixel_size):, int((window_size_y - CROP)/pixel_size):]
                print('wx', window_size_x, 'wy', window_size_y, 'CROP', crop, 'pixel_size', pixel_size)
                crop_point_x = int((crop)/pixel_size)
                crop_point_y = int((window_size_y - crop)/pixel_size)
                print(crop_point_x, crop_point_y)

                # stack = stack[:, :crop_point_x, crop_point_y:]
                stack = stack[:, int((window_size_x - crop)/pixel_size):, :int(crop/pixel_size)]
                # stack = stack[:, :int(CROP/pixel_size), :int(CROP/pixel_size)]
                # stack = stack[:, :int(CROP*pixel_size), :int(CROP*pixel_size)]

                # print()

            # stack = common.add_drift_intensity(stack, 1)

            print(stack.shape[1], 'x', stack.shape[2], 'px')


            stack = stack - stack.mean(axis=0) # remove space background

        if crop:
            in_crop = (particles[:, 0] < crop) & (particles[:, 1] < crop)
            particles = particles[in_crop, :]
            window_size_x = crop
            window_size_y = crop
        
        def add_outlines(timestep, ax):
            if particles.shape[1] == 4:
                particle_detection.show.add_particle_tracks(ax, particles, timestep, dimension=data_particles.get('dimension', 2),
                                                            window_size_x=window_size_x, window_size_y=window_size_y,)
            else:
                particle_detection.show.add_particle_outlines(
                    ax, particles, timestep, dimension=data_particles.get('dimension', 2),
                    outline=False, particle_diameter=data_particles['particle_diameter'],
                    window_size_x=window_size_x, window_size_y=window_size_y,)

        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, outfile,
                                                func=add_outlines, num_timesteps_in_data=num_timesteps,
                                                window_size_x=window_size_x, window_size_y=window_size_y,
                                                highlights=HIGHLIGHTS, every_nth_frame=every_nth_frame, output_type=output_type,
                                                annotation_color=annotation_color)
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