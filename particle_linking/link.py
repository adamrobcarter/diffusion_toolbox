import numpy as np
import pandas
import common
import trackpy
import tqdm

def go(file):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    pixel_size = data.get('pixel_size')
    particles  = data['particles']
    dimension = data.get('dimension', 2)

    assert particles.shape[1] == 3

    print('particles.dtype', particles.dtype)

    num_timesteps = np.unique(particles[:, dimension]).size
    print('av particles per frame', particles.shape[0]/num_timesteps)

    print('creating dataframe')

    dimension = data.get('dimension', 2)
    time_column = dimension
    id_column   = dimension + 1

    print('particles per frame', np.bincount(particles[:, time_column].astype('int')))


    if dimension == 2:
        columns = ['x', 'y', 'frame']
        columns_linked = ['x', 'y', 'frame', 'particle']
    else:
        columns = ['x', 'y', 'z', 'frame']
        columns_linked = ['x', 'y', 'z', 'frame', 'particle']

    if 'mixt' in file:
        # trackpy needs frame numbers not times
        # so we convert from time to number using the maps frame_to_time and time_to_frame
        frame_to_time = np.unique(particles[:, time_column])
        time_to_frame = lambda time: np.argmax(frame_to_time == time)

        for row in tqdm.trange(particles.shape[0], desc='mapping times', leave=False):
            particles[row, time_column] = time_to_frame(particles[row, time_column])

    # assert np.unique(particles[:, time_column]).size == particles[:, time_column].max() + 1, f'{np.unique(particles[:, time_column]).size} != {particles[:, time_column].max() + 1}'
    # if they are integers, this should check they're continuous and zero based

    features = pandas.DataFrame(particles, columns=columns)
    # the dtype of the frame column is float32
    print('created dataframe')
    loaded_df = False

    # print(features.dtypes)


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
    adaptive_stop = 1 * np.sqrt(2 * dt * D)
    memory = 1

    if file.startswith('eleanor'):
        pass
    if file.startswith('brennan'):
        pass
    if file.startswith('sim_nohydro'):
        pass
    if file.startswith('sim_hydro'):
        pass
    if file.startswith('sim_nointer'):
        pass
    if file == 'sophie1':
        search_range = 50
    if file.startswith('faxtor'):
        search_range = 15
    if file == 'carlos02':
        memory = 0

    print(f'search range = {search_range:.3g}um')
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

    print('trajectory lengths:')
    traj_lengths = np.bincount(trajs['particle'].to_numpy(dtype='int'))
    print(traj_lengths)
    common.term_hist(traj_lengths)

    print('filtering stubs')
    min_traj_length = 10
    # if file == 'sophie1':
    #     min_traj_length = 1
    num_trajs_before_filter = trajs.shape[0]
    trajs = trackpy.filter_stubs(trajs, min_traj_length)
    num_trajs = trajs.shape[0]
    print(f'dropped {num_trajs_before_filter-num_trajs} = {(num_trajs_before_filter-num_trajs)/num_trajs_before_filter*100:.0f}% of rows in filter_stubs')
    # filtering stubs might seem unneeded but it makes calculation of the MSD much much quicker

    print('trajectory lengths after filter:')
    traj_lengths = np.bincount(trajs['particle'].to_numpy(dtype='int'))
    print(traj_lengths)
    common.term_hist(traj_lengths, bins=np.arange(0, 100, 5))

    particles = trajs[columns_linked].to_numpy(dtype=particles.dtype)

    assert particles.shape[0] > 0, 'no particles in resulting dataset'

    if 'mixt' in file:
        # trackpy needs frame numbers not times
        # so we convert back from number to time
        for row in tqdm.trange(particles.shape[0], desc='remapping times', leave=False):
            particles[row, time_column] = frame_to_time[int(particles[row, time_column])]

    # after filtering stubs, the IDs are now no longer continuous
    IDs = np.unique(particles[:, id_column])
    print('num unique IDs', IDs.size)
    ID_map = {}
    for i, ID in enumerate(IDs):
        ID_map[ID] = i

        # # check that the time coordinate is continous while we're here
        # this_particle = particles[:, id_column] == ID
        # assert np.all(np.diff(np.unique(particles[this_particle, time_column])) == 1)

    for i in tqdm.trange(particles.shape[0], desc='updating IDs'):
        particles[i, id_column] = ID_map[particles[i, id_column]]

    
    if 'pack_frac_given' in data:
        exp_density = 4/np.pi * data['pack_frac_given'] / data['particle_diameter']**2
        exp_num = exp_density * data['window_size_x'] * data['window_size_y']
        print('num trajs', num_trajs, 'avg part per frame', particles.shape[0]/(particles[:, 2].max()+1), 'exp num', exp_num)

    
    num_timesteps = np.unique(particles[:, dimension]).size
    print('av particles per frame', particles.shape[0]/num_timesteps)

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

    # check continuity of time coordinate for each particle
    num_particles = np.unique(particles[:, id_column]).size
    bad = 0
    for ID in tqdm.trange(num_particles, desc='checking time continuity'):
        this_particle = particles[:, id_column] == ID
        times = np.unique(particles[this_particle, time_column])
        diffs = np.diff(times)
        if not np.all(diffs == 1):
            bad += 1
            print(times)

    assert bad / num_particles < 0.2, f'{bad} particles ({bad / num_particles:.1%}) had non-continuous time coordinate. Consider setting memory = 0'
    
    common.save_data(f'particle_linking/data/trajs_{file}.npz',
            particles=particles, radius=radius, time_step=data['time_step'],
            pixel_size=pixel_size,
            particle_diameter=data.get('particle_diameter'), pack_frac_given=data.get('pack_frac_given'),
            pack_frac=data.get('pack_frac'), particle_material=data.get('particle_material'),
            window_size_x=data.get('window_size_x'), window_size_y=data.get('window_size_y'),
            dimension=dimension)
            # particle_diameter_calced=particle_diameter_calced)

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file)