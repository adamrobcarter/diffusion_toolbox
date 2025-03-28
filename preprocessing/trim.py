import common

TRIMS = [1.0, 0.5, 0.25]
# TRIMS = [0.0625]

for file in common.files_from_argv('preprocessing/data', 'particles_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    for trim in TRIMS:
        newstack = stack[:int(stack.shape[0]*trim), :, :]
        newdata = dict(data)
        newdata['stack'] = newstack
        common.save_data(f'preprocessing/data/stack_{file}_trim{trim}.npz', **newdata)