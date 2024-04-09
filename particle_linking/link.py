import numpy as np
import pandas
import common
import trackpy

for file in common.files_from_argv('particle_detection/data', 'particlesdf_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    pixel_size = data['pixel_size']

    features = pandas.read_pickle(f'particle_detection/data/particlesdf_{file}.pickle')
    if file.startswith('eleanor'):
        search_range = 5
        filter_stubs = 1
    # search_range:
    # specify a maximum displacement, the farthest a particle can travel between frames. 
    # we should choose the smallest reasonable value because a large value slows computation time considerably.
    # filter_stubs:
    # ephemeral trajectories — seen only for a few frames — are usually spurious and never useful.
    # the convenience function filter_stubs keeps only trajectories that last for a given number of frames.
    trackpy.quiet()
    trajs = trackpy.link(features, search_range=search_range, memory=3) # memory lets particles dissapear for a frame
    print(trajs.describe())

    # num_trajs_before_filter = trajs['particle'].nunique()
    # trajs = trackpy.filter_stubs(trajs, )
    # num_trajs_after_filter = trajs['particle'].nunique()
    # print(f'dropped {(num_trajs_before_filter-num_trajs_after_filter)/num_trajs_before_filter:.2f} of trajs in filter')

    particles = trajs[['y', 'x', 'frame', 'particle']].to_numpy()
    particles[:, [0, 1]] *= pixel_size
    #                      ^    ^   I am aware these are the wrong way round
    # but it has to be so to work. possibly we introduced this error in the sparticles
    # tracking, but now we have to be consistant
    radius = trajs[['size']].to_numpy()[:, 0] # we use radius not diameter(size) for backward compatibility
    # ^^^ TODO I don't like this!!!!!!!!! call it diameter or size ffs
    
    np.savez(f'particle_linking/data/trajs_{file}.npz',
            #  particle_picklepath=picklepath,
            particles=particles, radius=radius, time_step=data['time_step'],
            pixel_size=pixel_size)
            # particle_diameter=particle_diameter,
            # particle_diameter_calced=particle_diameter_calced)