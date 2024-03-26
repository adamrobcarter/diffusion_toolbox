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

    print(stack.shape, stack.mean(axis=0).shape)
    stack = stack - stack.mean(axis=0) # remove space background
    print(stack.mean(axis=(1, 2)).std())
    # stack = stack - stack.mean(axis=(1, 2))[:, np.newaxis, np.newaxis] # remove total intensity fluctuations in time
        
    print(stack.mean(axis=(1, 2)).std())

    fig, ax = plt.subplots(1, 1)

    frames = range(0, min(stack.shape[0], 100), 1)

    def show(timestep):
        ax.clear()
        # plt.imshow(stack[timestep, :, :])
        # plt.imshow(stack[timestep, :, :])
        im = ax.imshow(stack[timestep, :, :], vmin=stack.min()/5, vmax=stack.max()/5)
        if timestep == 0:
            fig.colorbar(im)
        # plt.imshow(stack.min(axis=0))
        
        # print(stack[:, :, timestep].mean())

    common.save_gif(show, frames, fig, f"preprocessing/figures_png/stack_movie_{file}.gif", fps=8)