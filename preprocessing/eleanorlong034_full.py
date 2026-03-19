import numpy as np
import common
import os, time
import tqdm
import pandas

if __name__ == '__main__':
    total_rows = 0


    files = list(os.scandir('raw_data/eleanor/ref230815/patched/'))
    assert len(files)

    map = [np.array([])]*len(files)

    for f in tqdm.tqdm(files, desc='loading'):
        if f.name.startswith('2fps_s'):
            arr = pandas.read_pickle(f.path)

            number = int(f.name[6:-2]) - 1

            assert len(arr.index)
            total_rows += len(arr.index)

            map[number] = arr.to_numpy(dtype=np.float32)
            # map[number] = np.swapaxes(data, 1, 0)
        else:
            print('skipping', f.name)

    assert total_rows > 0

    print(total_rows, 'total_rows')
    particles = np.empty((total_rows, 4), dtype=np.float32)
    print(common.arraysize(particles))

    last_index = 0
    # last_time = 0

    for i, arr in enumerate(tqdm.tqdm(map, desc='combining data')):
        particles[last_index:last_index+arr.shape[0], :] = arr
        last_index += arr.shape[0]

    print(particles.shape)

    pixel_size = 0.288
    particles[:, [0, 1]] *= pixel_size

    num_timesteps = int(particles[:, 2].max()) + 1
    window_size_x = particles[:, 1].max()
    window_size_y = particles[:, 0].max()
    print(particles[:, 1].min(), particles[:, 1].max(), particles[:, 0].min(), particles[:, 0].max())

    EDGE_CROP = 4
    particles = common.crop_particles(particles, window_size_y-EDGE_CROP, window_size_x-EDGE_CROP, EDGE_CROP, EDGE_CROP)

    particle_diameter = 3.09

    metadata = dict(
        time_step = 0.5,
        pack_frac_given = 0.342,
        particle_diameter = particle_diameter,
    )

    common.save_trajs(
        'eleanorlong034',
        particles = particles,
        particles_labels = ['y', 'x', 't', 'id'],
        **metadata
    )
    common.save_particles(
        'eleanorlong034',
        particles = particles[:, [0, 1, 2]],
        particles_labels = ['y', 'x', 't'],
        **metadata
    )

    end_timestep = num_timesteps // 8
    particles = particles[particles[:, 2] < end_timestep, :]

    common.save_trajs(
        'eleanorlong034_div8',
        particles = particles,
        particles_labels = ['y', 'x', 't', 'id'],
        **metadata
    )
    common.save_particles(
        'eleanorlong034_div8',
        particles = particles[:, [0, 1, 2]],
        particles_labels = ['y', 'x', 't'],
        **metadata
    )



"""
The experiment reference (for me) is 230815.

This a zipped folder containing 145 files.
Compressed it’s 8GB and uncompressed it will be about 22GB.

These are python pickle files of Dataframes, managed with the pandas package.
e.g. data = pandas.read_pickle('filename.p')

Each file has four columns: 'y', 'x', 'frame', 'particle'.
 
'y' and 'x' give particle coordinates in pixels.
To convert to um use 0.288 (+/- 0.003) um/px
 
'frame' gives the time in frames.
To convert to seconds use 2fps (this is in all the file names by local convention).
 
'particle' gives an integer label of the particle identity.
 
The experiment is split into files each 1000 frames long.
Feel free to merge them if you have space in your RAM.
Particle identities are consistent from file to file.
    """