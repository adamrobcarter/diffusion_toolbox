import numpy as np
import matplotlib.pyplot as plt
import common
import sys

for file in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    # data = common.load('data/alice_stack_0.02.npz')
    # data = common.load('data/stack_alice_0.66.npz')
    stack = data['stack']

    print(stack.shape[1], 'x', stack.shape[2], 'px')

    for timestep in range(0, stack.shape[0], 4):
        plt.clf()
        # plt.imshow(stack[timestep, :, :])
        # plt.imshow(stack[timestep, :, :])

        
        # stack = stack - stack.mean(axis=0) # remove space background
        stack = stack.mean(axis=0) # remove space background

        # plt.imshow(stack[timestep, :, :])
        plt.imshow(stack[:, :])
        # plt.imshow(stack.min(axis=0))
        
        # excess = stack - stack.min(axis=0)
        # print(stack[:, :, timestep].mean())

        common.add_scale_bar(plt.gca(), data['pixel_size'], color='white')

        common.save_fig(plt.gcf(), f'preprocessing/figures_png/frame1_{file}.png', dpi=600, only_plot=True)
        break

        plt.pause(0.00001)