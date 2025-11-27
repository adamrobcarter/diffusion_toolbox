import common
import numpy as np
import matplotlib.pyplot as plt
import numba
import tqdm
import time

CROP = 0.2

# print('compiling')
# @numba.njit('float64[:,:](float64[:,:], float64[:,:], int32[:])', parallel=True)
@numba.njit(parallel=True)
def do_for_frame(this_frame, last_frame, x_min, x_max, y_min, y_max):
    corrs = np.full((x_max-x_min+1, y_max-y_min+1), np.nan)

    for ix, offset_x in enumerate(range(x_min, x_max+1)):
                
        if offset_x > 0:
            xstart_this = 0
            xend_this = -offset_x
            xstart_last = offset_x
            xend_last = last_frame.shape[0]
        elif offset_x < 0:
            xstart_this = -offset_x
            xend_this = this_frame.shape[0]
            xstart_last = 0
            xend_last = offset_x
        elif offset_x == 0:
            xstart_this = 0
            xend_this = this_frame.shape[0]
            xstart_last = 0
            xend_last = last_frame.shape[0]

        for iy, offset_y in enumerate(range(y_min, y_max+1)):
            if offset_y > 0:
                ystart_this = 0
                yend_this = -offset_y
                ystart_last = offset_y
                yend_last = last_frame.shape[1]
            elif offset_y < 0:
                ystart_this = -offset_y
                yend_this = this_frame.shape[1]
                ystart_last = 0
                yend_last = offset_y
            elif offset_y == 0:
                ystart_this = 0
                yend_this = this_frame.shape[1]
                ystart_last = 0
                yend_last = last_frame.shape[1]

            this_frame_moved   = this_frame[xstart_this:xend_this, ystart_this:yend_this]
            last_frame_cropped = last_frame[xstart_last:xend_last, ystart_last:yend_last]
                    
                # print(this_frame_moved.shape, last_frame_cropped.shape)
            corr = np.sum((this_frame_moved - last_frame_cropped)**2)
            corr /= this_frame_moved.size # compensate for the fact that we have different size arrays at different offsets
            # it would be better to crop for this surely!
            corrs[ix, iy] = corr

    assert np.isnan(corrs).sum() == 0
            
    return corrs
# print('done')

saved_offsets = {
    'psiche103': ((-8, 1), (-1, 32)),
    'psiche075': ((-4, 1), (-3, 1)),
}

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']
        pixel_size = data['pixel_size']

        # print('removing mean')
        # stack -= stack.mean(axis=0)

        border_x = int(stack.shape[1] * (1 - CROP) / 2)
        border_y = int(stack.shape[1] * (1 - CROP) / 2)

        s = stack.size
        # stack = stack[:, border_x:-border_x, border_y:-border_y]
        # stack = stack[:, ::8, ::8]
        print(f'smaller by {stack.size/s:.3f}')

        # frames_to_use = np.linspace(1, stack.shape[0]-1, 50)
        frames_to_use = range(stack.shape[0])

        all_corrs = np.full((len(frames_to_use), 2), np.nan)


        max_offset = 1

        # max_offsets = None
        # for key in saved_offsets.keys():
        #     if file.startswith(key):
        #         max_offsets = saved_offsets[key]
        #         print('got saved offsets')
        #         break
        # else:
        max_offsets = [[-max_offset, max_offset], [-max_offset, max_offset]]

        last = stack[0, :, :].astype(np.float32)

        fig, ax = plt.subplots(1, 1)

        for i, frame in enumerate(tqdm.tqdm(frames_to_use)):
            frame = int(frame)

            this = stack[frame, :, :].astype(np.float32)


            # print(this_frame, last_frame, max_offset)

            for repeat in range(50):
                corrs = do_for_frame(this, last, max_offsets[0][0], max_offsets[0][1], max_offsets[1][0], max_offsets[1][1])

                # if frame == 91:
                im = ax.imshow(corrs)
                    
                assert np.isnan(corrs).sum() == 0

                min = np.unravel_index(np.argmin(corrs, axis=None), corrs.shape) # https://numpy.org/doc/stable/reference/generated/numpy.argmax.html
                assert np.isnan(min).sum() == 0
                offset_px = np.array([max_offsets[0][0], max_offsets[1][0]]) + np.array(min)
                assert np.isnan(offset_px).sum() == 0
                all_corrs[i, :] = offset_px

                finished = True

                if max_offsets[0][0] >= offset_px[0]: # x min
                    max_offsets[0][0] = offset_px[0] - 1
                    print(f'min x updated to {max_offsets[0][0]}')
                    finished = False
                if max_offsets[1][0] >= offset_px[1]: # y min
                    max_offsets[1][0] = offset_px[1] - 1
                    print(f'min y updated to {max_offsets[1][0]}')
                    finished = False
                if max_offsets[0][1] <= offset_px[0]: # x max
                    max_offsets[0][1] = offset_px[0] + 1
                    print(f'max x updated to {max_offsets[0][1]}')
                    finished = False
                if max_offsets[1][1] <= offset_px[1]: # y max
                    max_offsets[1][1] = offset_px[1] + 1
                    print(f'max y updated to {max_offsets[1][1]}')
                    finished = False

                if finished:
                    break

            else:
                # raise Exception('I did not manage to find a minimum')
                break
            

        print(all_corrs[::20, :])
        assert np.isnan(all_corrs).sum() == 0
        print('x min max', all_corrs[:, 0].max(), all_corrs[:, 0].min())
        print('y min max', all_corrs[:, 1].max(), all_corrs[:, 1].min())
        common.save_data(f'preprocessing/data/corr_shift_{file}.npz', corrs=all_corrs)

        fig.colorbar(im)
        common.save_fig(fig, f'preprocessing/figures_png/correlations_{file}.png')

            # print(f'frame {frame} {(np.array(min))*pixel_size}um')
            # print(corrs.min(), corrs[min])
                    # print(corr)
            #         # print()

        