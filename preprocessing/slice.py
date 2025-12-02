import common
import numpy as np

SLICES = 5

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']

        slices = np.array_split(stack, SLICES, axis=1)
        [print(slice.shape) for slice in slices]

        for i, sliced_stack in enumerate(slices):
            common.save_data(
                f'preprocessing/data/stack_{file}_slice{i}',
                stack=sliced_stack,
                time_step=data['time_step'], pixel_size=data['pixel_size'],
                max_time_hours=data.get('max_time_hours'), particle_diameter=data.get('particle_diameter'),
                # window_size_(x/y) is not needed these days no?
            )