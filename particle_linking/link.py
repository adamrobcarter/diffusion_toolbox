import numpy as np
import pandas
import common
import trackpy

for file in common.files_from_argv('particle_detection/data', 'particlesdf_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    pixel_size    = data.get('pixel_size')
    window_size_x = data.get('window_size_x')
    window_size_y = data.get('window_size_y')


    try:
        features = pandas.read_pickle(f'particle_detection/data/particlesdf_{file}.pickle')
        loaded_df = True

    except FileNotFoundError:
        # we don't have the dataframe, so we gotta start from the raw particles
        particles = data['particles']
        features = pandas.DataFrame(particles, columns=['x', 'y', 'frame'])
        loaded_df = False

        # is there actually any point in the first method tbh?

    if file.startswith('eleanor'):
        search_range = 5
        filter_stubs = 1
        memory = 3
    if file.startswith('brennan'):
        search_range = 5 * 0.288
        filter_stubs = 1
        memory = 1
    # search_range:
    #   specify a maximum displacement, the farthest a particle can travel between frames. 
    #   we should choose the smallest reasonable value because a large value slows computation time considerably.
    # filter_stubs:
    #   ephemeral trajectories — seen only for a few frames — are usually spurious and never useful.
    #   the convenience function filter_stubs keeps only trajectories that last for a given number of frames.
    # memory:
    #   lets particles dissapear for a frame

    if features.shape[0] < 1e7:
        trackpy.quiet()
    trajs = trackpy.link(features, search_range=search_range, memory=memory)
    print(trajs.describe())

    # num_trajs_before_filter = trajs['particle'].nunique()
    # trajs = trackpy.filter_stubs(trajs, )
    # num_trajs_after_filter = trajs['particle'].nunique()
    # print(f'dropped {(num_trajs_before_filter-num_trajs_after_filter)/num_trajs_before_filter:.2f} of trajs in filter')

    particles = trajs[['y', 'x', 'frame', 'particle']].to_numpy()
    #                      ^    ^   I am aware these are the wrong way round
    # but it has to be so to work. possibly we introduced this error in the sparticles
    # tracking, but now we have to be consistant
    if pixel_size:
        particles[:, [0, 1]] *= pixel_size
    else:
        print('note: no pixel size')

    if loaded_df:
        radius = trajs[['size']].to_numpy()[:, 0] # we use radius not diameter(size) for backward compatibility
        # ^^^ TODO I don't like this!!!!!!!!! call it diameter or size ffs
    else:
        radius = None
    
    common.save_data(f'particle_linking/data/trajs_{file}.npz',
            #  particle_picklepath=picklepath,
            particles=particles, radius=radius, time_step=data['time_step'],
            pixel_size=pixel_size,
            particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'),
            num_timesteps=particles[:, 2].max(),
            window_size_x=window_size_x, window_size_y=window_size_y)
            # particle_diameter_calced=particle_diameter_calced)