import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

import preprocessing.stack_movie

USE_FIRST_FRAME_ONLY = False
FRAME_DIFF_TAU = 1

# REMOVE_BKG = True
# FRAME_DIFF = False

# FRAME_DIFF = True
# REMOVE_BKG = False

# FRAME_DIFF = False
# REMOVE_BKG = False

USE_ALL_FRAMES = True

def form_factor_sphere(q, R):
    Fs = 3 * (np.sin(q*R) - q*R*np.cos(q*R)) / (q*R)**3
    return Fs**2

def do_static_fourier(file, remove_bkg, frame_diff):
    
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']
    particle_diameter = data.get('particle_diameter')

    if remove_bkg:
        print('subtracting mean')
        stack = stack - stack.mean(axis=0)
        print(stack.shape, stack.mean(axis=0).shape)
        print(stack[0, :, :])
        # preprocessing.stack_movie.save_array_movie(stack, file=file, outputfilename='test.gif', pixel_size=pixel_size, time_step=1, window_size_x=stack.shape[1]*pixel_size, window_size_y=stack.shape[1]*pixel_size)

    print('min max mean', stack.min(), stack.max(), stack.mean())

    assert USE_ALL_FRAMES

    if USE_FIRST_FRAME_ONLY:
        if frame_diff:
            images = stack[[FRAME_DIFF_TAU], :, :] - stack[[0], :, :]
        else:
            images = stack[[0], :, :]
    else:
        if frame_diff:
            images = stack[FRAME_DIFF_TAU:, :, :] - stack[:-FRAME_DIFF_TAU, :, :]
        else:
            images = stack[:, :, :]
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
    title += f'_diff{FRAME_DIFF_TAU}' if frame_diff else ''
    title += '_bkgrem' if remove_bkg else ''

    common.save_data(f'DDM/data/static_fourier_{title}.npz',
                     fourier=fourier, fx=fx, fy=fy,
                     pixel_size=pixel_size, particle_diameter=particle_diameter,
                     num_frames_used=num_frames,
                     NAME=data.get('NAME'))

if __name__ == '__main__':
    errs = []

    for file in common.files_from_argv('preprocessing/data', 'stack_'):

        # try:
            do_static_fourier(file, False, False)
            do_static_fourier(file, True, False)
            do_static_fourier(file, False, True)
    #     except Exception as e:
    #         errs.append(e)

    # if len(errs):
    #     print(f'{len(errs)} exceptions')
    #     for e in errs:
    #         print(f"{type(e).__name__} at line {e.__traceback__.tb_lineno} of {__file__}: {e}")