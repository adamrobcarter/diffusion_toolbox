import common
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']


        jolt = 0.1

        newstack = np.zeros((stack.shape[0], int(stack.shape[1]*(1-jolt)), int(stack.shape[2]*(1-jolt))))

        newstack[:stack.shape[0]//2, :, :] = stack[:stack.shape[0]//2, :int(stack.shape[1]*(1-jolt)), :int(stack.shape[2]*(1-jolt))]

        newstack[stack.shape[0]//2:, :, :] = stack[stack.shape[0]//2:, int(stack.shape[1]*jolt):, :int(stack.shape[2]*(1-jolt))]

        common.save_data(f'preprocessing/data/stack_{file}_jolted.npz',
                        jolt=jolt,
                        pixel_size=data['pixel_size'], time_step=data['time_step'],
                        we_processed=data.get('we_processed'), nth_frame=data.get('nth_frame'))