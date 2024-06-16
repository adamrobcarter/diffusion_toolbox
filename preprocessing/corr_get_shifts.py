import common
import numpy as np
import matplotlib.pyplot as plt
import numba
import tqdm

# print('compiling')
# @numba.njit('float64[:,:](float64[:,:], float64[:,:], int32[:])', parallel=True)
@numba.njit(parallel=True)
def do_for_frame(this_frame, last_frame, max_offset):
    corrs = np.full((2*max_offset+1, 2*max_offset+1), np.nan)

    for ix, offset_x in enumerate(range(-max_offset, max_offset+1)):
        for iy, offset_y in enumerate(range(-max_offset, max_offset+1)):
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

            this_frame_moved   = this_frame[xstart_this:xend_this, ystart_this:yend_this]
            last_frame_cropped = last_frame[xstart_last:xend_last, ystart_last:yend_last]
                    
                # print(this_frame_moved.shape, last_frame_cropped.shape)
            corr = np.sum((this_frame_moved - last_frame_cropped)**2)
            corr /= this_frame_moved.size # compensate for the fact that we have different size arrays at different offsets
            corrs[ix, iy] = corr
    return corrs
# print('done')

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    pixel_size = data['pixel_size']

    # print('removing mean')
    # stack -= stack.mean(axis=0)

    # frames_to_use = np.linspace(1, stack.shape[0]-1, 50)
    frames_to_use = range(stack.shape[0])

    all_corrs = np.full((len(frames_to_use), 2), np.nan)


    max_offset = 35

    for i, frame in enumerate(tqdm.tqdm(frames_to_use)):
        frame = int(frame)

        this_frame = stack[frame,   :, :].astype(np.float64)
        last_frame = stack[0, :, :].astype(np.float64)
        this = np.zeros(this_frame.shape)
        this[:, :] = this_frame
        last = np.zeros(this_frame.shape)
        last[:, :] = last_frame
        # this_frame = stack[frame,   :, :]
        # last_frame = stack[0, :, :]


        # print(this_frame, last_frame, max_offset)
        corrs = do_for_frame(this, last, max_offset)

        # if frame == 91:
        #     im = plt.imshow(corrs)
        #     plt.colorbar(im)
        #     plt.show()
            
        assert np.isnan(corrs).sum() == 0

        min = np.unravel_index(np.argmin(corrs, axis=None), corrs.shape) # https://numpy.org/doc/stable/reference/generated/numpy.argmax.html
        offset_px = np.array(min) - max_offset
        all_corrs[i, :] = offset_px
        print(offset_px)
    
        assert  offset_px.max() < max_offset, 'you need a bigger max_offset'
        assert -offset_px.min() < max_offset, 'you need a bigger max_offset'

    print(all_corrs)
    print('max min', all_corrs.max(), all_corrs.min())
    common.save_data(f'preprocessing/data/corr_shift_{file}.npz', corrs=all_corrs)

        # print(f'frame {frame} {(np.array(min))*pixel_size}um')
        # print(corrs.min(), corrs[min])
                # print(corr)
        #         # print()

    