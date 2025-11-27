import common
import numpy as np
import scipy.ndimage
import skimage.filters
import tqdm

HIGHPASS_SIZE = 3 # px
# HIGHPASS_METHOD = 'gaussian'
HIGHPASS_METHOD = 'butterworth'

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']
        
        # stack = stack - stack.mean(axis=0)
        
        newdata = common.copy_not_stack(data)
        newdata['highpass_method'] = HIGHPASS_METHOD

        print('starting filter')

        if HIGHPASS_METHOD == 'gaussian':
            stack_lowpass = np.full_like(stack, np.nan)
            for frame in tqdm.trange(stack.shape[0]):
                stack_lowpass[frame] = scipy.ndimage.gaussian_filter(stack[frame], HIGHPASS_SIZE)
                newdata['highpass_size'] = HIGHPASS_SIZE

            stack_highpass = stack - stack_lowpass

        elif HIGHPASS_METHOD == 'butterworth':
            order = 10
            cutoff_px = 10
            cutoff = 1/cutoff_px * 0.5 # in [0, 0.5] where 0.5 is f_s/2
            stack_highpass = skimage.filters.butterworth(stack, cutoff_frequency_ratio=cutoff, high_pass=True, order=order, channel_axis=0).astype(stack.dtype)
            newdata['highpass_cutoff'] = cutoff
            newdata['highpass_order'] = order

        newdata['stack'] = stack_highpass
        print('finshed filter')

        common.save_data(f'preprocessing/data/stack_{file}_hpf.npz', **newdata)