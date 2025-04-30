import tifffile
import common
import sys

for file in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack'][:, :, :]
    stack -= stack.min()
    # print(stack.max())
    stack /= stack.max()
    stack *= 65535
    print(stack.shape)

    filename = f'preprocessing/data/stack_{file}.tiff'
    tifffile.imwrite(filename, stack.astype('uint16'))
    print(f'saved {filename}')