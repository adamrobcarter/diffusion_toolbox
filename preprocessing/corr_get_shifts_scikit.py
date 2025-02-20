import common
import numpy as np
import matplotlib.pyplot as plt
import numba
import tqdm
import time
import skimage

CROP = 0.2

saved_offsets = {
    'psiche103': ((-8, 1), (-1, 32)),
    'psiche075': ((-4, 1), (-3, 1)),
}

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    pixel_size = data['pixel_size']

    # print('removing mean')
    # stack -= stack.mean(axis=0)

    # border_x = int(stack.shape[1] * (1 - CROP) / 2)
    # border_y = int(stack.shape[1] * (1 - CROP) / 2)

    # s = stack.size
    # stack = stack[:, border_x:-border_x, border_y:-border_y]
    # stack = stack[:, ::8, ::8]
    # print(f'smaller by {stack.size/s:.3f}')

    # frames_to_use = np.linspace(1, stack.shape[0]-1, 50)
    frames_to_use = range(stack.shape[0])

    all_corrs = np.full((len(frames_to_use), 2), np.nan)


    # max_offset = 1

    # max_offsets = None
    # for key in saved_offsets.keys():
    #     if file.startswith(key):
    #         max_offsets = saved_offsets[key]
    #         print('got saved offsets')
    #         break
    # else:
    # max_offsets = [[-max_offset, max_offset], [-max_offset, max_offset]]

    last = stack[0, :, :].astype(np.float32)

    fig, ax = plt.subplots(1, 1)

    for i, frame in enumerate(tqdm.tqdm(frames_to_use)):
        frame = int(frame)

        this = stack[frame, :, :].astype(np.float32)

        shift, err, phasediff = skimage.registration.phase_cross_correlation(last, this)
        all_corrs[i, :] = shift
        print(shift)

        # print(this_frame, last_frame, max_offset)

        # for repeat in range(50):
        #     corrs = do_for_frame(this, last, max_offsets[0][0], max_offsets[0][1], max_offsets[1][0], max_offsets[1][1])

        #     # if frame == 91:
        #     im = ax.imshow(corrs)
                
        #     assert np.isnan(corrs).sum() == 0

        #     min = np.unravel_index(np.argmin(corrs, axis=None), corrs.shape) # https://numpy.org/doc/stable/reference/generated/numpy.argmax.html
        #     assert np.isnan(min).sum() == 0
        #     offset_px = np.array([max_offsets[0][0], max_offsets[1][0]]) + np.array(min)
        #     assert np.isnan(offset_px).sum() == 0
        #     all_corrs[i, :] = offset_px

        #     finished = True

        #     if max_offsets[0][0] >= offset_px[0]: # x min
        #         max_offsets[0][0] = offset_px[0] - 1
        #         print(f'min x updated to {max_offsets[0][0]}')
        #         finished = False
        #     if max_offsets[1][0] >= offset_px[1]: # y min
        #         max_offsets[1][0] = offset_px[1] - 1
        #         print(f'min y updated to {max_offsets[1][0]}')
        #         finished = False
        #     if max_offsets[0][1] <= offset_px[0]: # x max
        #         max_offsets[0][1] = offset_px[0] + 1
        #         print(f'max x updated to {max_offsets[0][1]}')
        #         finished = False
        #     if max_offsets[1][1] <= offset_px[1]: # y max
        #         max_offsets[1][1] = offset_px[1] + 1
        #         print(f'max y updated to {max_offsets[1][1]}')
        #         finished = False

        #     if finished:
        #         break

        # else:
        #     # raise Exception('I did not manage to find a minimum')
        #     break
        

    print(all_corrs[::20, :])
    assert np.isnan(all_corrs).sum() == 0
    print('x min max', all_corrs[:, 0].max(), all_corrs[:, 0].min())
    print('y min max', all_corrs[:, 1].max(), all_corrs[:, 1].min())
    common.save_data(f'preprocessing/data/corr_shift_scikit_{file}.npz', corrs=all_corrs)

    # fig.colorbar(im)
    # common.save_fig(fig, f'preprocessing/figures_png/correlations_{file}.png')

        # print(f'frame {frame} {(np.array(min))*pixel_size}um')
        # print(corrs.min(), corrs[min])
                # print(corr)
        #         # print()

    