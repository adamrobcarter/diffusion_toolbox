import numpy as np
import pandas
import common
import trackpy
import tqdm

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    pixel_size = data.get('pixel_size')
    particles  = data['particles']

    assert particles.shape[1] == 3

    print('particles.dtype', particles.dtype)

    print('creating dataframe')
    features = pandas.DataFrame(particles, columns=['x', 'y', 'frame'])
    # the dtype of the frame column is float32
    print('created dataframe')
    loaded_df = False

    print(features.dtypes)


    # steps are distributed ~ N(0, sqrt(2 D dt)^2)
    # so search_range = sqrt(2 D dt) would get 68% of steps
    # 2 * sqrt(2 D dt) : 95%
    # 3 * sqrt(2 D dt) : 99.7%
    # however for some reason at the mo I seem to need more like 10 *, idk why
    dt = data['time_step']
    if 'eleanor' in file or 'sim_nohydro' in file or 'sim_hydro' in file or 'brennan' in file:
        D = 0.04
        if data['pack_frac_given'] > 0.3:
            D = 0.03
        if data['pack_frac_given'] > 0.6:
            D = 0.02
    else:
        D = 1
    search_range = 10 * np.sqrt(2 * dt * D)
    adaptive_stop = 2 * np.sqrt(2 * dt * D)
    memory = 1

    if file.startswith('eleanor'):
        pass
    if file.startswith('brennan'):
        pass
    if file.startswith('sim_nohydro'):
        pass
    if file.startswith('sim_nointer'):
        pass
    print('search range', search_range)
    # search_range:
    #   specify a maximum displacement, the farthest a particle can travel between frames. 
    #   we should choose the smallest reasonable value because a large value slows computation time considerably.
    # filter_stubs:
    #   ephemeral trajectories — seen only for a few frames — are usually spurious and never useful.
    #   the convenience function filter_stubs keeps only trajectories that last for a given number of frames.
    # memory:
    #   lets particles dissapear for a frame

    # if features.shape[0] < 1e7:
    trackpy.quiet()

    print('beginning trackpy')
    trajs = trackpy.link(features, search_range=search_range, memory=memory, adaptive_stop=adaptive_stop)
    print(trajs.describe())

    print('filtering stubs')
    num_trajs_before_filter = trajs.shape[0]
    trajs = trackpy.filter_stubs(trajs, 10)
    num_trajs = trajs.shape[0]
    print(f'dropped {(num_trajs_before_filter-num_trajs)/num_trajs_before_filter:.2f} of rows in filter_stubs')
    # filtering stubs might seem unneeded but it makes calculation of the MSD much much quicker

    particles = trajs[['y', 'x', 'frame', 'particle']].to_numpy(dtype=particles.dtype)
    #                   ^    ^   I am aware these are the wrong way round
    # but it has to be so to work. possibly we introduced this error in the sparticles
    # tracking, but now we have to be consistant

    # after filtering stubs, the IDs are now no longer continuous
    IDs = np.unique(particles[:, 3])
    ID_map = {}
    for i, ID in enumerate(IDs):
        ID_map[ID] = i
    for i in tqdm.trange(particles.shape[0], desc='updating IDs'):
        particles[i, 3] = ID_map[particles[i, 3]]
    
    if 'pack_frac_given' in data:
        exp_density = 4/np.pi * data['pack_frac_given'] / data['particle_diameter']**2
        exp_num = exp_density * data['window_size_x'] * data['window_size_y']
        print('num trajs', num_trajs, 'avg part per frame', particles.shape[0]/(particles[:, 2].max()+1), 'exp num', exp_num)

    
    num_before = particles.shape[0]

    print('particles dtype', particles.dtype)

    # num_after = particles.shape[0]
    # print(f'rows before:{num_before}, rows after:{num_after}')
    # assert num_after > 0.8*num_before

    # check trajectory lengths
    # traj_lengths = np.full((num_trajs,), np.nan)
    # for i in tqdm.trange(num_trajs, desc='getting trajectory lengths'):
    #     traj_lengths[i] = np.sum(particles[:, 3] == i)
    # print('traj lengths:')
    # common.term_hist(traj_lengths)

    if loaded_df:
        radius = trajs[['size']].to_numpy(dtype=particles.dtype)[:, 0] # we use radius not diameter(size) for backward compatibility
        # ^^^ TODO I don't like this!!!!!!!!! call it diameter or size ffs
    else:
        radius = None
    
    common.save_data(f'particle_linking/data/trajs_{file}.npz',
            particles=particles, radius=radius, time_step=data['time_step'],
            pixel_size=pixel_size,
            particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'),
            window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'))
            # particle_diameter_calced=particle_diameter_calced)