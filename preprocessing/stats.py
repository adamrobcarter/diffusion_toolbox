# this file is just for seeing what's saved inside a data file

import common

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')