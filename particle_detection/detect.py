import common
import numpy as np
import stracking.detectors
import time
import sys

for datasource in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{datasource}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    particle_diameter = data['particle_diameter']
    time_step         = data['time_step']
    num_timesteps = stack.shape[0]

    frame_average = stack.mean(axis=0)
    stack = stack - frame_average[np.newaxis, :, :] # could subtract the pixel min not avg (or median)

    overlap = 0.5
    if datasource == 'pierre_simdownsampled':
        min_sigma = 0.5
        max_sigma = 2
        threshold = 0.1
    if datasource == 'pierre_sim':
        min_sigma = 2
        max_sigma = 8
        threshold = 0.1
    elif datasource == 'pierre_exp': 
        min_sigma = 3
        max_sigma = 5
        threshold = 3
    elif datasource == 'eleanor0.34':
        min_sigma = 2
        max_sigma = 8
        threshold = 25
    elif datasource == 'eleanor0.01':
        min_sigma = 2
        max_sigma = 8
        threshold = 25
    elif datasource.startswith('marine'):
        min_sigma = 1
        max_sigma = 4
        threshold = 1500
        overlap = 0.5

    t0 = time.time()

    print('starting detector')

    detector = stracking.detectors.DoGDetector(min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    detector_output = detector.run(stack)

    print(f'detector finished in {time.time()-t0:.0f}s')

    particles = detector_output.data # rows of t, x, y
    # if datasource == 'sim':
    #     PIXEL = 0.1
    # if datasource == 'sim_downsampled':
    #     PIXEL = 0.4
    particles[:, [1, 2]] *= pixel_size

    # rearrange from t, x, y to x, y, t
    particles = particles[:, [1, 2, 0]]

    print(detector_output.properties.keys())
    # print(detector_output.properties['radius'].shape)
    # print(detector_output.properties['radius'].mean(), detector_output.properties['radius'].std(), detector_output.properties['radius'].max())

    # picklepath = f'data/removed_av_particles_{datasource}.pickle'
    # with open(picklepath, 'wb') as handle:
    #     pickle.dump(detector_output, handle)

    radius = detector_output.properties['radius']

    particle_diameter_calced = 2 * np.sqrt(2) * radius.mean() * pixel_size
        # there is a line in the DoGDetector source about this sqrt 2
    print(f'calced diameter {particle_diameter_calced:.3f}um')

    # print(f'found {(particles.shape[0]/num_timesteps):0f} particles per frame')
    # density = particles.shape[0]/num_timesteps / (stack.shape[0]*stack.shape[1]*pixel_size**2)
    # pack_frac = np.pi/4 * density * particle_diameter**2
    # print(f'so packing fraction phi = {pack_frac:.3f}')

    np.savez(f'particle_detection/data/particles_{datasource}.npz',
            #  particle_picklepath=picklepath,
            particles=particles, radius=radius, time_step=time_step,
            min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold, computation_time=time.time()-t0,
            pixel_size=pixel_size, num_timesteps=num_timesteps, particle_diameter=particle_diameter,
            particle_diameter_calced=particle_diameter_calced)