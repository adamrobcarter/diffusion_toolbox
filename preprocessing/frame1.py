import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm

for file in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    # data = common.load('data/alice_stack_0.02.npz')
    # data = common.load('data/stack_alice_0.66.npz')
    stack = data['stack']

    print(stack.shape[1], 'x', stack.shape[2], 'px')


        
    # stack = stack - stack.mean(axis=0) # remove space background
    # stack = stack.mean(axis=0)
    plt.figure(figsize=(3, 3))

    plt.imshow(stack[0, :, :], cmap=matplotlib.cm.Greys)
    # plt.imshow(stack.min(axis=0))
    
    # excess = stack - stack.min(axis=0)
    # print(stack[:, :, timestep].mean())

    plt.ylim(000, 600)
    plt.xlim(000, 600)

    common.add_scale_bar(plt.gca(), data['pixel_size'], color='white')

    common.save_fig(plt.gcf(), f'preprocessing/figures_png/frame1_{file}.png', dpi=600, only_plot=True)
    common.save_fig(plt.gcf(), f'/home/acarter/presentations/intcha24/figures/frame1_{file}.png', dpi=600, only_plot=True)