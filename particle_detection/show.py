import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys
import warnings
import preprocessing.stack_movie


def add_particle_outlines(ax, pixel_size, particles, radius, timestep, channel=None, outline=True):
    # radius can be None

    particles_at_t = particles[:, 2] == timestep
    # print(particles_at_t.sum())
    if particles_at_t.sum() == 0:
        warnings.warn(f'No particles found at timestep')
    # assert particles_at_t.sum() > 0
    # plt.scatter(particles[particles_at_t, X_INDEX]/pixel_size, particles[particles_at_t, Y_INDEX]/pixel_size, s=50*radius[particles_at_t]**2*pixel_size**2,
    #             facecolors='none', edgecolors='red', alpha=0.5, linewidth=0.8)

    # x = particles[particles_at_t, 1]/pixel_size
    # y = particles[particles_at_t, 0]/pixel_size
    x = particles[particles_at_t, 1]
    y = particles[particles_at_t, 0]

    # r = radius[particles_at_t] * np.sqrt(2) # TODO: should this be /pixel_size?
    r = np.full_like(x, 40*pixel_size**2)
    # warnings.warn('i disabled showing radius')
    if particles.shape[1] == 4:
        id = particles[particles_at_t, 3]

    def color(i):
        if particles.shape[1] == 4:
            return matplotlib.cm.tab20(int(id[i]%20))
        else:
            if channel == 'red':
                return 'white'
            elif channel == 'green':
                return 'red'
            return 'red'
        
    alpha = 0.7

    # cross = False
    # if cross:
    #     pass
    if outline:
        circles = [plt.Circle((x[i], y[i]), edgecolor=color(i), facecolor='none', radius=r[i]*2, linewidth=1, alpha=alpha) for i in range(particles_at_t.sum())]
        # c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)
        c = matplotlib.collections.PatchCollection(circles, match_original=True)
    else:
        circles = [plt.Circle((x[i], y[i]), facecolor=color(i), radius=r[i], alpha=alpha) for i in range(particles_at_t.sum())]
        c = matplotlib.collections.PatchCollection(circles, match_original=True)

    ax.add_collection(c)

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):

        data_particles = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data_particles['particles']
        radius    = data_particles.get('radius', np.zeros(particles.shape[0]))
        time_step = data_particles['time_step']
        channel   = data_particles.get('channel')

        try:
            data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
            stack      = data_stack['stack']
            pixel_size = data_stack['pixel_size']
            
            # crop
            # stack = stack[:, :500, :500]

            stack = stack - stack.mean(axis=0) # remove space background
            # stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

        except FileNotFoundError:
            num_timesteps = int(particles[:, 2].max() - 1)
            stack = np.zeros((num_timesteps, 320, 320))
            pixel_size = 1
            # radius = np.full(particles.shape[0], 0.002*160)
            radius = np.full(particles.shape[0], 0.002*160*pixel_size**2/10)
            print()

        def add_outlines(timestep, ax):
            add_particle_outlines(ax, pixel_size, particles, radius, timestep, channel=channel, outline=False)

        figsize = 0.01 * np.array([data_particles['window_size_x'], data_particles['window_size_y']])

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_ylim(particles[:, 0].min(), particles[:, 0].max())
        ax.set_xlim(particles[:, 1].min(), particles[:, 1].max())
        ax.set_aspect('equal')

        add_outlines(0, ax)

        common.save_fig(fig, f'particle_detection/figures_png/frame1_{file}.png')