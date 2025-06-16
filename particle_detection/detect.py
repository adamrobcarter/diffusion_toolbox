import common
import numpy as np
import time
import sys
import trackpy
import matplotlib.pyplot as plt
import termplotlib
import tqdm

NO_MINMASS = False

if NO_MINMASS:
    print('warning, running with no minmass')

for file in common.files_from_argv('preprocessing/data/', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    particle_diameter = data.get('particle_diameter')
    time_step         = data.get('time_step')
    depth_of_field    = data.get('depth_of_field')
    num_timesteps = stack.shape[0]

    if num_timesteps > 1:
        stack = stack - stack.mean(axis=0)
    print('stack dtype', stack.dtype)

    print('converting stack range')
    stack = np.interp(stack, (stack.min(), stack.max()), (0, 1))
    # convert to 0->1 range. why are we doing this?
    print(stack.max(), stack.mean(), 'maxmean', stack.min(), stack.std())

    print('stack dtype', stack.dtype)

    threshold = 1/255
    maxsize = None
    minmass = 0
    diameter = None
    separation = None # min separation between features. default diameter + 1
    invert = False
    percentile = None # "Features must have a peak brighter than pixels in this percentile. This helps eliminate spurious peaks. default 64

    if file == 'pierre_simdownsampled':
        pass
    elif file == 'pierre_sim':
        pass
    elif file == 'pierre_exp': 
        diameter = 7
        minmass = 0.2
        separation = 5
        # would threshold be good here, as it's meant to be for when the bkg is noisy?
    elif file == 'eleanor0.34' or file == 'eleanor0.01':
        diameter = 7
        minmass = 1
        
    elif file == 'eleanor0.34' or file == 'eleanor0.01' or file == 'eleanor0.34_ds1' or file == 'eleanor0.01_ds1':
        diameter = 7
        minmass = 1
        
    elif file == 'eleanor0.34_ds2' or file == 'eleanor0.01_ds2':
        diameter = 5
        minmass = 1/2
        
    elif file == 'eleanor0.34_ds4' or file == 'eleanor0.01_ds4':
        diameter = 3
        minmass = 0.01
        separation = 2.8
        
    elif file == 'eleanor0.34_ds8' or file == 'eleanor0.01_ds8':
        diameter = 3
        minmass = 0.0001
        separation = 2.8
        percentile = 20

    elif file.startswith('marine2'):
        diameter = 5
        minmass = 0.01

    elif file.startswith('marine3'):
        diameter = 9
        minmass = 0.0001
        percentile = 20

    elif file.startswith('marine'):
        # for marine10 at least
        diameter = 5
        minmass = 0.2

    elif file == 'victor0':
        diameter = 15
        minmass = 2.8
        invert = True
        threshold *= 10

    elif file == 'carlos01':
        diameter = 7
        # minmass = 2.8
        # invert = Tru
        # threshold *= 10

    else:
        raise Exception('you need to provide parameters for this dataset')

    if NO_MINMASS:
        minmass = 0

    outfile = f'{file}_nominmass' if NO_MINMASS else file

    t0 = time.time()

    print('starting detector')

    assert diameter is not None, "you haven't yet provided `diameter`"

    # features = trackpy.batch(stack[:100], 11, processes=16)
    trackpy.quiet() # turns off reporting every frame as it's processed
    progress = tqdm.tqdm(total=stack.shape[0])
    def update(frame, features):
        progress.update()
        return features
    features = trackpy.batch(stack, processes=16, diameter=diameter, minmass=minmass,
                             separation=separation, threshold=threshold, invert=invert,
                             after_locate=update)
    progress.close()

    print(features.describe())

    print('mass:')
    counts, bin_edges = np.histogram(features['mass'], bins=20)
    term_fig = termplotlib.figure()
    term_fig.hist(counts, bin_edges, force_ascii=False, orientation="horizontal")
    term_fig.show()

    # hist_fig, hist_ax = plt.subplots(1, 1, figsize=(3, 3))
    # hist_ax.hist(features['mass'], bins=np.linspace(0, np.max(features['mass']), 20))
    # hist_ax.set_yticks([])
    # hist_ax.set_xlabel('mass')
    # hist_ax.set_ylabel('count')
    # common.save_fig(hist_fig, f'/home/acarter/presentations/cin_first/figures/mass_hist_{outfile}.png', hide_metadata=True)
    # common.save_fig(hist_fig, f'particle_detection/figures_png/mass_hist_{outfile}.png')

    # detector = stracking.detectors.DoGDetector(min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    # detector_output = detector.run(stack)

    print(f'detector finished in {time.time()-t0:.0f}s')

    particles = features[['y', 'x', 'frame']].to_numpy()
    #                      ^    ^   I am aware these are the wrong way round
    # but it has to be so to work. possibly we introduced this error in the sparticles
    # tracking, but now we have to be consistant

    there is still a switch here which is annoying
    is this connected to needing the right orientation for faxtor?
    it would be really good to properly fix this

    particles[:, [0, 1]] *= pixel_size

    radius = features[['size']].to_numpy()[:, 0] # we use radius not diameter(size) for backward compatibility

    # plt.hist(features[['mass']].to_numpy(), bins=20)
    # plt.semilogy()
    # plt.savefig('hist.png')

    particle_diameter_calced = 2 * radius.mean() * pixel_size
        # there is a line in the DoGDetector source about this sqrt 2
    print(f'calced diameter {particle_diameter_calced:.3f}um')

    assert particles.shape[0] > 0, 'no particles were found'

    to_save = dict(
        #  particle_picklepath=picklepath,
        particles=particles, radius=radius, time_step=time_step,
        threshold=threshold, maxsize=maxsize, minmass=minmass, diameter=diameter, separation=separation,
        computation_time=time.time()-t0, depth_of_field=depth_of_field,
        pixel_size=pixel_size, num_timesteps=num_timesteps, particle_diameter=particle_diameter,
        particle_diameter_calced=particle_diameter_calced,
        window_size_x=pixel_size*stack.shape[1], window_size_y=pixel_size*stack.shape[2],
        channel=data.get('channel'), NAME=data.get('NAME'),
        pack_frac_given=data.get('pack_frac_given'),
    )

    print(f'found {(particles.shape[0]/num_timesteps):.0f} particles per frame')
    if particle_diameter is not None and not np.isnan(particle_diameter):
        pack_frac = common.calc_pack_frac(particles, particle_diameter, stack.shape[0]*pixel_size, stack.shape[1]*pixel_size)
        print(f'so packing fraction phi = {pack_frac:.3f}')

        to_save['pack_frac'] = pack_frac

    common.save_data(f'particle_detection/data/particles_{outfile}.npz', **to_save)