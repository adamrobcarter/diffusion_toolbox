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
    stack = stack - stack.mean(axis=0)
        
    fig, ax = plt.subplots(1, 1)

    frames = range(0, min(stack.shape[0], 100), 1)

    def show(timestep):
        ax.clear()
        # plt.imshow(stack[timestep, :, :])
        # plt.imshow(stack[timestep, :, :])
        ax.imshow(stack[timestep, :, :])
        # plt.imshow(stack.min(axis=0))
        
        # print(stack[:, :, timestep].mean())

    common.save_gif(show, frames, fig, f"preprocessing/figures_png/stack_movie_{file}.gif", fps=8)