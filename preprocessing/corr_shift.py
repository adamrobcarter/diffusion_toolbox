import common
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']

        for frame in range(1, stack.shape[0]):

            this_frame = stack[frame,   :, :]
            last_frame = stack[frame-1, :, :]

            max_offset = 10

            corrs = np.full((2*max_offset+1, 2*max_offset+1), np.nan)

            for ix, offset_x in enumerate(range(-max_offset, max_offset)):
                for iy, offset_y in enumerate(range(-max_offset, max_offset)):
                    if offset_x > 0:
                        if offset_y > 0:
                            this_frame_moved   = this_frame[offset_x:, offset_y:]
                            last_frame_cropped = last_frame[:offset_x, :offset_y]
                        else:
                            this_frame_moved   = this_frame[offset_x:, :offset_y]
                            last_frame_cropped = last_frame[:offset_x, offset_y:]
                    else:
                        if offset_y > 0:
                            this_frame_moved   = this_frame[:offset_x, offset_y:]
                            last_frame_cropped = last_frame[offset_x:, :offset_y]
                        else:
                            this_frame_moved   = this_frame[:offset_x, :offset_y]
                            last_frame_cropped = last_frame[offset_x:, offset_y:]

                    corr = np.sum(this_frame_moved * last_frame_cropped)
                    corrs[ix, iy] = corr
                    print(corrs)