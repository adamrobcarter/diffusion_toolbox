import common
import numpy as np
import matplotlib.pyplot as plt
import numba
import tqdm

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):

        data2 = common.load(f'preprocessing/data/corr_shift_{file}.npz')
        corrs = data2['corrs']
        
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']
        pixel_size = data['pixel_size']

        lost_x = int(corrs[:, 0].max() - corrs[:, 0].min())
        lost_y = int(corrs[:, 1].max() - corrs[:, 1].min())
        
        # strategy is to line them all up in an array padded by nan, then crop
        newstack = np.full((stack.shape[0], stack.shape[1]+lost_x, stack.shape[2]+lost_y), np.nan, dtype=np.float16)

        for t in tqdm.trange(stack.shape[0]):
            x_offset = int(corrs[t, 0] - corrs[:, 0].min())
            y_offset = int(corrs[t, 1] - corrs[:, 1].min())
            newstack[t, x_offset:x_offset+stack.shape[1], y_offset:y_offset+stack.shape[2]] = stack[t, :, :]

        # now crop
        newstack = newstack[:, lost_x:-lost_x, lost_y:-lost_y]

        assert not np.any(np.isnan(newstack))
        assert newstack.size

        common.save_data(f'preprocessing/data/stack_{file}_shifted.npz',
                    stack=newstack,
                    pixel_size=data['pixel_size'], time_step=data['time_step'],
                    we_processed=data.get('we_processed'), nth_frame=data.get('nth_frame'))


