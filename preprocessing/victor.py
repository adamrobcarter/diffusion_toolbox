import imageio
import numpy as np
import common

for file in ['victor0', 'victor0_scalebar']:
    image = imageio.imread(f'raw_data/victor/{file}.png')

    # this image is already greyscale so we just take the first channel
    image = image[:, :, 0]

    # give it a reduntant time dimension
    stack = image[np.newaxis, :, :,]

    common.save_data(f'preprocessing/data/stack_{file}.npz', stack=stack, pixel_size=0.0005)