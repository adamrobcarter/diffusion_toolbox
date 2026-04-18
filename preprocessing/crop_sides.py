import common

CROP = 256

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']

        s = stack.shape[1]

        newstack = stack[:, CROP:-CROP, :]

        newdata = common.copy_not_stack(data)
        newdata['stack'] = newstack

        common.save_data(f'preprocessing/data/stack_{file}_cropsides{CROP}.npz', **newdata)
