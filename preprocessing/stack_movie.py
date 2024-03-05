import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.animation
import tqdm

for file in sys.argv[1:]:
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    # data = common.load('data/alice_stack_0.02.npz')
    # data = common.load('data/stack_alice_0.66.npz')
    stack = data['stack']

    print(stack.shape[1], 'x', stack.shape[2], 'px')

        
    fig, ax = plt.subplots(1, 1)

    frames = range(0, stack.shape[0], 1)
    progress = tqdm.tqdm(total=len(frames))

    def show(timestep):
        ax.clear()
        # plt.imshow(stack[timestep, :, :])
        # plt.imshow(stack[timestep, :, :])
        ax.imshow(stack[timestep, :, :]-stack.min(axis=0))
        # plt.imshow(stack.min(axis=0))
        
        # excess = stack - stack.min(axis=0)
        # print(stack[:, :, timestep].mean())
        progress.update()

    ani = matplotlib.animation.FuncAnimation(fig, show, frames=frames, interval=500, repeat=False)
    ani.save(f"preprocessing/figures_png/stack_movie_{file}.gif", dpi=300,
            writer=matplotlib.animation.PillowWriter(fps=8))
    progress.close()