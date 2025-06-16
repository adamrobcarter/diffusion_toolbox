import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft
import warnings
import DDM.static_fourier_show

import preprocessing.stack_movie

USE_FIRST_FRAME_ONLY = False
FRAME_DIFF_TAU = 1

# REMOVE_BKG = True
# FRAME_DIFF = False

FRAME_DIFF = True
REMOVE_BKG = False

# FRAME_DIFF = False
# REMOVE_BKG = False

USE_ALL_FRAMES = True

NORMALISE_INTENSITY_IN_TIME = False

def form_factor_sphere(q, R):
    Fs = 3 * (np.sin(q*R) - q*R*np.cos(q*R)) / (q*R)**3
    return Fs**2

def do_static_fourier(file, remove_bkg, frame_diff):
    
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    # stack      = data['stack'].astype(np.float32)
    stack      = data['stack']
    pixel_size = data['pixel_size']
    particle_diameter = data.get('particle_diameter')

    assert np.isfinite(stack).all()

    if remove_bkg:
        print('subtracting mean')
        stack = stack - stack.mean(axis=0)
        # print(stack.shape, stack.mean(axis=0).shape)
        # print(stack[0, :, :])
        # preprocessing.stack_movie.save_array_movie(stack, file=file, outputfilename='test.gif', pixel_size=pixel_size, time_step=1, window_size_x=stack.shape[1]*pixel_size, window_size_y=stack.shape[1]*pixel_size)

    # print('min max mean', stack.min(), stack.max(), stack.mean())

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

    if frame_diff and images.dtype != np.float32:
        warnings.warn('should use float32 when doing diff?')

    assert np.isfinite(images).all()

    if NORMALISE_INTENSITY_IN_TIME:
        avg_intensity_per_frame = images.mean(axis=(1, 2))
        # print(avg_intensity_per_frame)
        # assert np.all(avg_intensity_per_frame != 0)
        # assert np.isfinite(avg_intensity_per_frame).all()
        images = images / avg_intensity_per_frame[:, np.newaxis, np.newaxis]

        bad = ~ np.isfinite(images)
        print('bad', bad.sum(), bad.sum()/bad.size)

        images[bad] = 1

        assert np.isfinite(images).all()
    
    # image = image[::4, ::4]
    # pixel_size *= 4

    # rng = np.random.default_rng()
    # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
    #                     2 * np.pi * np.arange(200) / 34)
    # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

    assert np.isfinite(images).all()

    if False:
        assert images.mean() != 0
        images = images / images.mean()
        
    assert np.isfinite(images).all()

    fx, fy, fouriers = common.fourier_2D(images, pixel_size, (1, 2))
    del images
    # print(fouriers.shape)
    # fourier = fouriers.mean(axis=0)
    # fx = fx.mean(axis=0)
    # fy = fy.mean(axis=0)
    # print('f', fourier.shape, fx.shape)
    # print('fourier mean', (np.abs(fourier)**2).mean())

    fouriers_mod_sq = np.abs(fouriers)**2
    fourier_mod_sq = fouriers_mod_sq.mean(axis=0)

    title = file
    title += f'_diff{FRAME_DIFF_TAU}' if frame_diff else ''
    title += '_bkgrem' if remove_bkg else ''

    common.save_data(f'DDM/data/static_fourier_{title}.npz',
                     fourier_mod_sq=fourier_mod_sq, fx=fx, fy=fy,
                     pixel_size=pixel_size, particle_diameter=particle_diameter,
                     num_frames_used=num_frames, time_step=data['time_step'],
                     NAME=data.get('NAME'))

if __name__ == '__main__':
    errs = []

    for file in common.files_from_argv('preprocessing/data', 'stack_'):

        # try:
            do_static_fourier(file, remove_bkg=REMOVE_BKG, frame_diff=FRAME_DIFF)

            title = file
            title += f'_diff{FRAME_DIFF_TAU}' if FRAME_DIFF else ''
            title += '_bkgrem' if REMOVE_BKG else ''
            DDM.static_fourier_show.go(title)

            # do_static_fourier(file, remove_bkg=False, frame_diff=True)
            # do_static_fourier(file, remove_bkg=True,  frame_diff=False)
            # do_static_fourier(file, remove_bkg=False, frame_diff=False)
    #     except Exception as e:
    #         errs.append(e)

    # if len(errs):
    #     print(f'{len(errs)} exceptions')
    #     for e in errs:
    #         print(f"{type(e).__name__} at line {e.__traceback__.tb_lineno} of {__file__}: {e}")