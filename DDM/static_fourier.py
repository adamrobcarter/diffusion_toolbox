import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

USE_FIRST_FRAME_ONLY = False
# FRAME_DIFF = True
# REMOVE_BKG = False
FRAME_DIFF = False
REMOVE_BKG = True
USE_ALL_FRAMES = True

def form_factor_sphere(q, R):
    Fs = 3 * (np.sin(q*R) - q*R*np.cos(q*R)) / (q*R)**3
    return Fs**2

def do_static_fourier(file):
    
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']
    particle_diameter = data.get('particle_diameter')

    if REMOVE_BKG:
        print('subtracting mean')
        stack = stack - stack.mean(axis=0)

    assert USE_ALL_FRAMES

    if USE_FIRST_FRAME_ONLY:
        if FRAME_DIFF:
            images = stack[[1], :, :] - stack[[0], :, :]
        else:
            images = stack[[0], :, :]
    else:
        if FRAME_DIFF:
            images = stack[1::, :, :] - stack[:-1:, :, :]
        else:
            images = stack[::, :, :]
    del stack
    num_frames = images.shape[0]
    
    # image = image[::4, ::4]
    # pixel_size *= 4

    # rng = np.random.default_rng()
    # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
    #                     2 * np.pi * np.arange(200) / 34)
    # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

    print(images.shape)
    fx, fy, fouriers = common.fourier_2D(images, pixel_size, (1, 2))
    del images
    # print(fouriers.shape)
    fourier = fouriers.mean(axis=0)
    # fx = fx.mean(axis=0)
    # fy = fy.mean(axis=0)
    print('f', fourier.shape, fx.shape)

    

    title = file
    title += '_diff' if FRAME_DIFF else ''
    title += '_bkgrem' if REMOVE_BKG else ''

    common.save_data(f'DDM/data/static_fourier_{title}.npz',
                     fourier=fourier, fx=fx, fy=fy,
                     pixel_size=pixel_size, particle_diameter=particle_diameter,
                     num_frames_used=num_frames,
                     NAME=data.get('NAME'))

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):

        do_static_fourier(file)