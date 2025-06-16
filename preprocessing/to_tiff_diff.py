import tifffile
import common
import sys

for file in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack'][:, :, :]

    diffs = stack[1:, :, :] - stack[:-1, :, :]
    # stack -= stack.min()
    # # print(stack.max())
    # stack /= stack.max()
    # stack *= 65535
    
    filename = f'preprocessing/tiffs/stack_{file}_diff.tiff'
    tifffile.imwrite(filename, diffs)#.astype('uint16'))
    print(f'saved {filename}')