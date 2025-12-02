import common
import scipy.ndimage
import numpy as np

LENGTH = 50

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')

        stack = data['stack']

        print('finding background')
        bkg = scipy.ndimage.uniform_filter(stack.astype(np.float32), size=(LENGTH, 0, 0), mode='nearest').astype(np.float16)
        # note in the docs there is a warning about inaccuracies if you use a limited-precision data type
        print('todo: if still using float32, keep it in float32 uptil here?')
        
        print('subtracting background')
        stack -= bkg

        print('saving')
        newdata = common.copy_not_stack(data)
        newdata['stack'] = stack[LENGTH//2:-LENGTH//2, :, :] # remove start and end cause we can't filter properly there
        newdata['moving_avg_length'] = LENGTH

        common.save_data(f'preprocessing/data/stack_{file}_movavrem.npz', **newdata)