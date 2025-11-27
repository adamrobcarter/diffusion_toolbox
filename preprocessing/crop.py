import common

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']

        s = stack.shape[1]

        newstack = stack[:, s//3:2*s//3, s//3:2*s//3]

        newdata = dict(data)
        newdata['stack'] = newstack

        common.save_data(f'preprocessing/data/stack_{file}_crop.npz', **newdata)