import common

LENGTH = 100
START = 0

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')

        stack = data['stack']

        newdata = common.copy_not_stack(data)
        newdata['stack'] = stack[START:START+LENGTH, :, :]

        common.save_data(f'preprocessing/data/stack_{file}_extract{START}_{LENGTH}.npz', **newdata)