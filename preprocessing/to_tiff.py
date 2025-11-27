import tifffile
import common
import sys

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack'][:, :, :]
        # stack -= stack.min()
        # # print(stack.max())
        # stack /= stack.max()
        # stack *= 65535
        print(stack.shape)

        filename = f'preprocessing/tiffs/stack_{file}.tiff'
        tifffile.imwrite(filename, stack)#.astype('uint16'))
        print(f'saved {filename}')