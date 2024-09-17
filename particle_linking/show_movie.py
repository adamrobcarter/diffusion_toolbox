import common
import preprocessing.stack_movie
import particle_detection.show_movie
import sys

for file in sys.argv[1:]:
    try:
        data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data_stack['stack']
        pixel_size = data_stack['pixel_size']
        time_step  = data_stack['time_step']
        found_stack = True
    except FileNotFoundError:
        found_stack = False

    data_particles = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data_particles['particles']
    radius    = data_particles['radius']

    # crop
    stack = stack[:, :500, :500]

    stack = stack - stack.mean(axis=0) # remove space background
    
    def add_outlines(timestep, ax):
        particle_detection.show.add_particle_outlines(ax, pixel_size, particles, radius, timestep)

    filename = f'movie_linked_{file}'

    if found_stack:
        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"particle_linking/figures_png/{filename}.gif",
                                                func=add_outlines)
    else:
        preprocessing.stack_movie.save_array_movie(time_step=time_step, file=file, outputfilename=f"particle_linking/figures_png/{filename}.gif",
                                                func=add_outlines, no_stack=True, num_frames=50)
    # preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4",
    #                                            func=add_outlines)