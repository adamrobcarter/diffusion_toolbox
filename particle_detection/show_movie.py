import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.animation
import tqdm
import particle_detection.show


def go(filename, stack_file, particles_file, output_file):
    # this is also used by particle_linking.show_movie
    
    data = common.load(stack_file)
    stack             = data['stack']
    pixel_size        = data['pixel_size']
    # particle_diameter = data['particle_diameter']
    # num_timesteps = stack.shape[0]

    # print(stack.shape[1], 'x', stack.shape[2], 'px')

    stack = stack - stack.mean(axis=0)
    stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

    # print('stack min max', stack.min(), stack.max())

    data = common.load(particles_file)
    particles = data['particles']
    radius    = data['radius']
    
    print(f'found {(particles.shape[0]/particles[:, 2].max()):0f} particles per frame')
    # print(particles[:, 2].min(), particles[:, 2].max())

    fig, ax = plt.subplots(1, 1)

    end_frame = min(stack.shape[0], 50)
    frames = range(0, end_frame, 1)

    def show(timestep):
        ax.clear()

        particle_detection.show.show_frame(fig, ax, stack, pixel_size, particles, radius, timestep, filename)
        ax.set_title(f'{filename} t={timestep:03d}')

    common.save_gif(show, frames, fig, output_file, fps=2)

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        go(file,
           stack_file=f'preprocessing/data/stack_{file}.npz',
           particles_file=f'particle_detection/data/particles_{file}.npz',
           output_file=f"particle_detection/figures_png/movie_{file}.gif")