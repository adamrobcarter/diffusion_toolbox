"""
Extract a subset of frames from the given stack and save as a new stack.

Usage:
    ``python -m preprocessing.extract --start=START --length=LENGTH dataset1 [dataset2 ...]``

Options:
    --start=START       Starting frame index
    --length=LENGTH     Number of frames to extract
"""

import common

if __name__ == '__main__':

    argparser = common.argparser()
    argparser.add_argument('--start',  type=int, help='Starting frame index')
    argparser.add_argument('--length', type=int, help='Number of frames to extract')
    args = argparser.parse_args()

    for file in args.files:
        data = common.load(f'preprocessing/data/stack_{file}.npz')

        stack = data['stack']

        newdata = common.copy_not_stack(data)
        newdata['stack'] = stack[args.start:args.start+args.length, :, :]

        common.save_data(f'preprocessing/data/stack_{file}_extract{args.start}_{args.length}.npz', **newdata)