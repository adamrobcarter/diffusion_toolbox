import czifile
import common
import numpy as np

for i, file in enumerate(common.get_directory_files('raw_data/marine', 'czi')):
    image = czifile.imread(file)
    stack = image.squeeze()
    print(f'{file} is marine{i}', stack.shape)
    np.savez(f'preprocessing/data/stack_marine{i}', stack=stack)