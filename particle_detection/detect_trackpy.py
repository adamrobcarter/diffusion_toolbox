import common
import numpy as np
import time
import sys
import trackpy
import matplotlib.pyplot as plt

for datasource in common.files_from_argv('preprocessing/data/', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{datasource}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    particle_diameter = data['particle_diameter']
    time_step         = data['time_step']
    depth_of_field    = data.get('depth_of_field')
    num_timesteps = stack.shape[0]

    stack = stack - stack.mean(axis=0)
    stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

    threshold = 1/255
    maxsize = None
    minmass = 0
    diameter = None
    separation = None

    if datasource == 'pierre_simdownsampled':
        pass
        # min_sigma = 0.5
        # max_sigma = 2
        # threshold = 0.1
    if datasource == 'pierre_sim':
        pass
        # min_sigma = 2
        # max_sigma = 8
        # threshold = 0.1
    elif datasource == 'pierre_exp': 
        diameter = 7
        minmass = 0.2
        separation = 5
        # would threshold be good here, as it's meant to be for when the bkg is noisy?
    elif datasource == 'eleanor0.34' or datasource == 'eleanor0.01':
        diameter = 7
        minmass = 2
        # min_sigma = 2
        # max_sigma = 8
        # threshold = 25

    elif datasource.startswith('marine'):
        # for marine10 at least
        diameter = 5
        minmass = 0.2
        # min_sigma = 1
        # max_sigma = 4
        # threshold = 1500
        # overlap = 0.5

    t0 = time.time()

    print('starting detector')

    assert diameter is not None, "you haven't yet provided `diameter`"

    # features = trackpy.batch(stack[:100], 11, processes=16)
    trackpy.quiet() # turns off reporting every frame as it's processed
    features = trackpy.batch(stack, processes=16, diameter=diameter, minmass=minmass,
                             separation=separation, threshold=threshold)
    print(features.describe())

    # detector = stracking.detectors.DoGDetector(min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    # detector_output = detector.run(stack)

    print(f'detector finished in {time.time()-t0:.0f}s')

    particles = features[['y', 'x', 'frame']].to_numpy()
    #                      ^    ^   I am aware these are the wrong way round
    # but it has to be so to work. possibly we introduced this error in the sparticles
    # tracking, but now we have to be consistant

    # radius = 
    print(particles.shape)
    # if datasource == 'sim':
    #     PIXEL = 0.1
    # if datasource == 'sim_downsampled':
    #     PIXEL = 0.4
    particles[:, [0, 1]] *= pixel_size

    # print(detector_output.properties['radius'].shape)
    # print(detector_output.properties['radius'].mean(), detector_output.properties['radius'].std(), detector_output.properties['radius'].max())

    # picklepath = f'data/removed_av_particles_{datasource}.pickle'
    # with open(picklepath, 'wb') as handle:
    #     pickle.dump(detector_output, handle)

    radius = features[['size']].to_numpy()[:, 0] # we use radius not diameter(size) for backward compatibility
    print('radius shape', radius.shape)

    plt.hist(features[['mass']].to_numpy(), bins=20)
    plt.semilogy()
    plt.savefig('hist.png')

    particle_diameter_calced = 2 * radius.mean() * pixel_size
        # there is a line in the DoGDetector source about this sqrt 2
    print(f'calced diameter {particle_diameter_calced:.3f}um')

    print(f'found {(particles.shape[0]/num_timesteps):0f} particles per frame')
    density = particles.shape[0]/num_timesteps / (stack.shape[0]*stack.shape[1]*pixel_size**2)
    pack_frac = np.pi/4 * density * particle_diameter**2
    print(f'so packing fraction phi = {pack_frac:.3f}')
    
    print('it', particles[:, 2].min(), particles[:, 2].max())

    # np.savez(f'particle_detection/data/particles_{datasource}_trackpy.npz',
    np.savez(f'particle_detection/data/particles_{datasource}.npz',
            #  particle_picklepath=picklepath,
            particles=particles, radius=radius, time_step=time_step,
            threshold=threshold, maxsize=maxsize, minmass=minmass, diameter=diameter, separation=separation,
            computation_time=time.time()-t0, depth_of_field=depth_of_field,
            pixel_size=pixel_size, num_timesteps=num_timesteps, particle_diameter=particle_diameter,
            particle_diameter_calced=particle_diameter_calced, method='trackpy')