import tifffile
import common
import sys

for file in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack'][1, :, :]
    stack -= stack.min()
    print(stack.max())
    stack *= 65535

    tifffile.imwrite(f'preprocessing/data/stack_{file}.tiff', stack.astype('uint16'))