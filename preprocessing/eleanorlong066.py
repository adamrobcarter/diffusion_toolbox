import numpy as np
import common
import os, time
import tqdm

if __name__ == '__main__':
    total_rows = 0

    map = [np.array([])]*155

    for f in tqdm.tqdm(list(os.scandir('raw_data/eleanor/eleanorlong066')), desc='loading'):
        if f.name.startswith('s_'):
            data = np.load(f.path)

            number = int(f.name[2:-4]) - 2
            
            total_rows += data.shape[1]
            map[number] = np.swapaxes(data, 1, 0)

    print(total_rows, 'total_rows')
    particles = np.empty((total_rows, 3), dtype=np.float32)
    print(common.arraysize(particles))

    last_index = 0
    # last_time = 0

    for i, arr in enumerate(tqdm.tqdm(map, desc='combining data')):
        data = arr.astype(np.float32)
        # first_time = data[0, 2]
        # print('times', last_time, first_time, data[-1, 2])

        data[:, 2] += i*1000
        particles[last_index:last_index+arr.shape[0], :] = data
        last_index += arr.shape[0]

        # last_time = data[-1, 2]

    print(particles.shape)

    particles[:, 2] -= particles[:, 2].min() # make zero-based

    pixel_size = 0.288
    particles[:, [0, 1]] *= pixel_size

    num_timesteps = int(particles[:, 2].max()) + 1
    window_size_x = particles[:, 0].max()
    window_size_y = particles[:, 1].max()

    EDGE_CROP = 4
    particles = common.crop_particles(particles, window_size_x-EDGE_CROP, window_size_y-EDGE_CROP, EDGE_CROP, EDGE_CROP)
    window_size_x -= 2*EDGE_CROP
    window_size_y -= 2*EDGE_CROP

    particle_diameter = 3.09

    t0 = time.time()
    np.save(f'particle_detection/data/particles_eleanorlong066.npy', particles)
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')

    t0 = time.time()
    common.save_data(f'particle_detection/data/particles_eleanorlong066.npz', particles=particles,
            time_step=0.5, particle_diameter=particle_diameter, pixel_size=pixel_size,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.656)
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')

    end_timestep = num_timesteps // 2
    particles = particles[particles[:, 2] < end_timestep, :]

    t0 = time.time()
    common.save_data(f'particle_detection/data/particles_eleanorlong066_div2.npz', particles=particles,
            time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.656)
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')

    end_timestep = num_timesteps // 4
    particles = particles[particles[:, 2] < end_timestep, :]

    t0 = time.time()
    common.save_data(f'particle_detection/data/particles_eleanorlong066_div4.npz', particles=particles,
            time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.656)
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')

    end_timestep = num_timesteps // 8
    particles = particles[particles[:, 2] < end_timestep, :]

    t0 = time.time()
    common.save_data(f'particle_detection/data/particles_eleanorlong066_div8.npz', particles=particles,
            time_step=0.5, particle_diameter=2.82, pixel_size=pixel_size,
            window_size_x=window_size_x, window_size_y=window_size_y,
            pack_frac_given=0.656)
    t1 = time.time()
    print(f'took {t1-t0:.0f}s')



    """
    This should contain x,y,t for the duration of the high phi experiment. It is 
    saved in 155 separate numpy arrays (columns x,y,t) each spanning 1000 frames 
    of the experiment (this is an artefact of the way we save videos but I have 
    found helpful for dealing with larger data sets piecewise rather than all at 
    once). The slightly annoying thing is that to cut it down to 16 bit precision 
    the t coordinate ranges from 0-999 in each file, but it should be okay to add 
    the right multiple of 1000 based on the index in the file name.

    I don't know how much ram you have - but if you can't get it into a format 
    where all the coordinates are in a single file I can send you the hacky code 
    I use for counting in boxes one file at a time and then stitching it all 
    together at the end

    Sorry this is slightly hacky! Hope it's not too irritating to use. 21hrs of 
    high packing fraction is just a lot of coordinates haha
    """