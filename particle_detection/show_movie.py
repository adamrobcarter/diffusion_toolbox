import numpy as np
import common
import sys
import particle_detection.show
import preprocessing.stack_movie
import typing

CROP = None
HIGHLIGHTS = True
FORCE_FPS = 10

def go(file, infile=None, outfile=None, crop=False, every_nth_frame=None,
       output_type : typing.Literal['movie', 'frames'] = 'movie',
       dpi=300, tracks=False, highlights=False, show_blobs_too=False,
       particle_color='red',
        **kwargs):
        """
        output_type: 'movie' will produce a gif, 'frames' saves each frame of the movie to a separate .png
        tracks: True will plot the history of the trajectory not just the current position

        annotation_color should be in kwargs if needed
        """
        
        if not infile:
            infile = f'particle_detection/data/particles_{file}.npz'
        assert outfile
        
        if 'trajs' in infile and tracks == False:
            print('WARNING the colours seem to be messed up when tracks=False')

        data_particles = common.load(infile)
        particles = data_particles['particles']
        # particles = particles[:, [1, 0, 2]] # DON'T ASK ME WHY

        print('particles max', particles[:, 0].max(), particles[:, 0].min())
        print('particles max', particles[:, 1].max(), particles[:, 1].min())
        
        time_step = data_particles['time_step']
        time_column = data_particles.get('dimension', 2)
        num_timesteps = particles[:, time_column].max() + 1
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
            if tracks:
                assert particles is not None
                particle_detection.show.add_particle_tracks(ax, particles, timestep, dimension=data_particles.get('dimension', 2),
                                                            window_size_x=window_size_x, window_size_y=window_size_y)
            if not tracks or show_blobs_too:
                particle_detection.show.add_particle_outlines(
                    ax, particles, timestep, dimension=data_particles.get('dimension', 2),
                    outline=False, particle_diameter=data_particles.get('particle_diameter', 5),
                    window_size_x=window_size_x, window_size_y=window_size_y, color=particle_color)

        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, outfile,
                                                func=add_outlines, num_timesteps_in_data=num_timesteps,
                                                window_size_x=window_size_x, window_size_y=window_size_y,
                                                highlights=highlights, every_nth_frame=every_nth_frame, output_type=output_type,
                                                dpi=dpi, force_fps=FORCE_FPS,
                                                **kwargs)


if __name__ == '__main__':

    for file in sys.argv[1:]:
        highlights = '_highlights' if HIGHLIGHTS else ''
        go(
            file,
            outfile = f'particle_detection/figures_png/movie_{file}{highlights}.gif',
            highlights=HIGHLIGHTS
        )