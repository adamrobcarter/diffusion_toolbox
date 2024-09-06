import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys
import warnings
import preprocessing.stack_movie


def add_particle_outlines(ax, pixel_size, particles, radius, timestep, channel=None):
    # radius can be None

    particles_at_t = particles[:, 2] == timestep
    # print(particles_at_t.sum())
    if particles_at_t.sum() == 0:
        warnings.warn(f'No particles found at timestep')
    # assert particles_at_t.sum() > 0
    # plt.scatter(particles[particles_at_t, X_INDEX]/pixel_size, particles[particles_at_t, Y_INDEX]/pixel_size, s=50*radius[particles_at_t]**2*pixel_size**2,
    #             facecolors='none', edgecolors='red', alpha=0.5, linewidth=0.8)

    x = particles[particles_at_t, 1]/pixel_size
    # common.term_hist(x)
    y = particles[particles_at_t, 0]/pixel_size
    r = radius[particles_at_t] * np.sqrt(2) # TODO: should this be /pixel_size?
    r = np.full_like(r, r.mean())
    warnings.warn('i disabled showing radius')
    if particles.shape[1] == 4:
        id = particles[particles_at_t, 3]

    def color(i):
        if particles.shape[1] == 4:
            return matplotlib.cm.tab20(int(id[i]%20))
        else:
            if channel:
                return 'white'
            return 'red'
        
    alpha = 0.7

    cross = False
    outline = True
    if cross:
        pass
    elif outline:
        circles = [plt.Circle((x[i], y[i]), edgecolor=color(i), facecolor='none', radius=r[i]*2, linewidth=0.5, alpha=alpha) for i in range(particles_at_t.sum())]
        # c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)
        c = matplotlib.collections.PatchCollection(circles, match_original=True)
    else:
        circles = [plt.Circle((x[i], y[i]), facecolor=color(i), radius=r[i], alpha=alpha) for i in range(particles_at_t.sum())]
        c = matplotlib.collections.PatchCollection(circles, match_original=True)

    ax.add_collection(c)

if __name__ == '__main__':
    for file in sys.argv[1:]:
        data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data_stack['stack']
        pixel_size = data_stack['pixel_size']
        time_step  = data_stack['time_step']
        channel    = data_stack.get('channel')

        data_particles = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data_particles['particles']
        radius    = data_particles['radius']

        # crop
        # stack = stack[:, :500, :500]

        stack = stack - stack.mean(axis=0) # remove space background
        
        def add_outlines(timestep, ax):
            add_particle_outlines(ax, pixel_size, particles, radius, timestep, channel=channel)

        filename = f'movie_{file}'
        preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"particle_detection/figures_png/{filename}.gif",
                                                func=add_outlines, channel=channel)
        # preprocessing.stack_movie.save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4",
        #                                            func=add_outlines)

"""
if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack             = data['stack']
        pixel_size        = data['pixel_size']
        particle_diameter = data.get('particle_diameter')
        num_timesteps = stack.shape[0]
        
        if num_timesteps > 1:
            stack = stack - stack.mean(axis=0)
        stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

        # frame_average = stack.mean(axis=0)
        # stack = stack - frame_average[np.newaxis, :, :]
        # stack += frame_average / 10
        # stack = frame_average / 10 - stack
        
        # stack = stack - stack.mean(axis=0) # remove space background
        
        # stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        radius = data['radius']
        print(particles.shape)

        ##### radii are wrong! pls fix!!

        if particle_diameter is not None:
            density = particles.shape[0]/num_timesteps / (stack.shape[0]*stack.shape[1]*pixel_size**2)
            pack_frac = np.pi/4 * density * particle_diameter**2
            print(f'packing fraction phi = {pack_frac:.3f}')

        print(f'radius min', radius.min(), 'median', np.median(radius), 'max', radius.max())
        print(f'radius min', radius.min()* np.sqrt(2), 'median', np.median(radius)* np.sqrt(2), 'max', radius.max()* np.sqrt(2))

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        # fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        fig.tight_layout()

        # for timestep in range(0, 1000, 10):
        if True:
            timestep = num_timesteps // 2

            show_frame(fig, ax, stack, pixel_size, particles, radius, timestep, file)

            common.save_fig(fig, f'particle_detection/figures_png/particles_{file}.png', dpi=300)

            """