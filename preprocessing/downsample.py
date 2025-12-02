import common
import skimage.measure
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']
        print(stack.shape)
        width = stack.shape[1]
        height = stack.shape[2]

        for downsampling in [1, 2, 4, 8]:
            new_width = width // downsampling
            new_height = height // downsampling

            trimmed_stack = stack[:, :new_width*downsampling, :new_height*downsampling] # this removes a few rows/columns such that the width and height are multiples of downsampling

            downsampled_stack = skimage.measure.block_reduce(trimmed_stack, block_size=(1, downsampling, downsampling), func=np.mean)
            downsampled_stack = downsampled_stack.astype(np.uint8)

            common.save_data(f'preprocessing/data/stack_{file}_ds{downsampling}.npz',
                            stack=downsampled_stack, time_step=data['time_step'], pixel_size=data['pixel_size']*downsampling,
                            pack_frac_given=data.get('pack_frac_given'), particle_diameter=data.get('particle_diameter'))