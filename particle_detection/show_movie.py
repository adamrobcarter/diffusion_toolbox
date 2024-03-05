import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.animation
import tqdm
import particle_detection.show

for file in sys.argv[1:]:
    
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    particle_diameter = data['particle_diameter']
    num_timesteps = stack.shape[0]

    print(stack.shape[1], 'x', stack.shape[2], 'px')

    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    radius = data['radius']
    print(particles.shape)

    fig, ax = plt.subplots(1, 1)

    end_frame = min(stack.shape[0], 100)
    frames = range(0, end_frame, 1)
    # progress = tqdm.tqdm(total=len(frames))

    def show(timestep):
        ax.clear()

        particle_detection.show.show_frame(fig, ax, stack, pixel_size, particles, radius, timestep)
        ax.set_title(f'{file} t={timestep:03d}')

        # progress.update()

    # ani = matplotlib.animation.FuncAnimation(fig, show, frames=frames, interval=500, repeat=False)
    # ani.save(f"particle_detection/figures_png/movie_{file}.gif", dpi=300,
    #         writer=matplotlib.animation.PillowWriter(fps=2))
    # progress.close()
    common.save_gif(show, frames, fig, f"particle_detection/figures_png/movie_{file}.gif", fps=2)