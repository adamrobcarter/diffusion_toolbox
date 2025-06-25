import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import tqdm
import preprocessing.stack_movie
import scattering_functions

# CROP = 150
# X_START = 250
# Y_START = 50
CROP = False


def find_color(row, bhwindow=False, window_size_x=None, window_size_y=None, channel=None):
    if row.size == 4:
        c = matplotlib.cm.tab20(int(row[3]%20))
    else:
        if channel == 'red':
            c = 'grey' # was white
        elif channel == 'green':
            c = 'red'
        else:
            c = 'red'

    if bhwindow:
        alpha = scattering_functions.blackman_harris_window(window_size_x, window_size_y, row[0], row[1])
    else:
        alpha = 1

    return c, alpha

def add_particle_outlines(ax, particles, timestep, dimension, particle_diameter, channel=None, outline=True,
                          window_size_x=None, window_size_y=None, bhwindow=False):
    # radius can be None

    assert particles.size

    particles_at_t = particles[:, dimension] == timestep
    # print(particles_at_t.sum())
    if particles_at_t.sum() == 0:
        raise Exception(f'No particles found at timestep {timestep}')
    # assert particles_at_t.sum() > 0
    # plt.scatter(particles[particles_at_t, X_INDEX]/pixel_size, particles[particles_at_t, Y_INDEX]/pixel_size, s=50*radius[particles_at_t]**2*pixel_size**2,
    #             facecolors='none', edgecolors='red', alpha=0.5, linewidth=0.8)

    x = window_size_x - particles[particles_at_t, 0]
    # i think this probably happens because when we plot the image the origin is bottom left
    # but when trackpy does the positions it's top left, or something
    y = particles[particles_at_t, 1]

    radius = particle_diameter / 2

    alpha_mult = 1
    color, alpha = find_color(particles[particles_at_t, :][0, :], bhwindow=bhwindow, window_size_x=window_size_x, window_size_y=window_size_y, channel=channel)

    # cross = False
    # if cross:
    #     pass
    if outline:
        circles = [plt.Circle((x[i], y[i]), edgecolor=color, facecolor='none', radius=radius*2, linewidth=3, alpha=alpha*alpha_mult) for i in range(particles_at_t.sum())]
        # c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)
        assert len(circles)
        c = matplotlib.collections.PatchCollection(circles, match_original=True)
    else:
        circles = [plt.Circle((x[i], y[i]), facecolor=color, radius=radius, alpha=alpha*alpha_mult) for i in range(particles_at_t.sum())]
        assert len(circles)
        c = matplotlib.collections.PatchCollection(circles, match_original=True)

    ax.add_collection(c)

def add_particle_tracks(ax, particles, timestep, dimension, window_size_x, window_size_y):
    assert particles.size

    # before_t = particles[:, dimension] <= timestep
    
    # if before_t.sum() == 0:
    #     raise Exception(f'No particles found at timestep {timestep}')
    
    for i, id in enumerate(np.unique(particles[:, dimension + 1])):
    # for i, id in enumerate(tqdm.tqdm(np.unique(particles[:, dimension + 1]))):

        particles_this_id = particles[particles[:, 3] == id, :]
        particles_this_time = particles_this_id[particles_this_id[:, 2] <= timestep, :]
        x = window_size_x - particles_this_time[:, 0]# - X_START
        ## ^^^^^^^^^^^^^^ this again si because trackpy uses top left as origin
        # not bottom left, really you should correct that in the particle detection
        y = particles_this_time[:, 1]# - Y_START

        ax.plot(x, y, linewidth=2)
        # ax.plot(x, 1024 - y) # please don't ask me

        # if i > 500:
        #     break
    
    # # particles_before

    # x = particles[particles_at_t, 0]
    # y = particles[particles_at_t, 1]

    # # r = radius[particles_at_t] * np.sqrt(2) # TODO: should this be /pixel_size?
    # r = np.full_like(x, 1.5) # you would lose this problem if u actually showed the radius u numpty
    # print('U HARDCODED THE ABOVE NUMBER U NUMPTY')
    # # warnings.warn('i disabled showing radius')
    # if dimension == 2 and particles.shape[1] > 3:
    #     id = particles[particles_at_t, 3]
    # if dimension == 3 and particles.shape[1] > 4:
    #     id = particles[particles_at_t, 4]

    # radius = particle_diameter / 2

    # alpha_mult = 1

    # # cross = False
    # # if cross:
    # #     pass
    # if outline:
    #     circles = [plt.Circle((x[i], y[i]), edgecolor=color(i)[0], facecolor='none', radius=radius*2, linewidth=3, alpha=color(i)[1]*alpha_mult) for i in range(particles_at_t.sum())]
    #     # c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)
    #     assert len(circles)
    #     c = matplotlib.collections.PatchCollection(circles, match_original=True)
    # else:
    #     circles = [plt.Circle((x[i], y[i]), facecolor=color(i)[0], radius=radius, alpha=color(i)[1]*alpha_mult) for i in range(particles_at_t.sum())]
    #     assert len(circles)
    #     c = matplotlib.collections.PatchCollection(circles, match_original=True)

    # ax.add_collection(c)

def go(file, ax, fig, bhwindow=False):
    # need to pass `fig` so that we can resize it to fit the size of the data
    
    filename_particles = file
    if file.endswith('_first'):
        filename_particles = file.split('_first')[0]

    data_particles = common.load(f'particle_detection/data/particles_{filename_particles}.npz')
    particles = data_particles['particles']
    radius    = data_particles.get('radius', np.zeros(particles.shape[0]))
    time_step = data_particles['time_step']
    channel   = data_particles.get('channel')

    try:
        data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data_stack['stack']
        pixel_size = data_stack['pixel_size']
            
            # crop
            # stack = stack[:, X_START:X_START+CROP, Y_START:Y_START+CROP]
            # stack = stack[:, :500, :500]
        print('removing background')
        stack = stack - stack.mean(axis=0) # remove space background
            # stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range
        no_stack = False

    except FileNotFoundError:
        no_stack = True
        num_timesteps = int(particles[:, 2].max() - 1)
            # stack = np.zeros((num_timesteps, 320, 320)) # this is very silly
        pixel_size = 1
        print()

        # particles = particles[:, [1, 0, 2]] # DON'T ASK ME WHY


    def add_outlines(timestep, ax):
        add_particle_outlines(ax, particles, timestep, dimension=data_particles.get('dimension', 2), particle_diameter=data_particles.get('particle_diameter', 3),
                                  channel=channel, outline=False, window_size_x=data_particles['window_size_x'], window_size_y=data_particles['window_size_y'],
                                  bhwindow=bhwindow)

    figsize = 4 * np.array([data_particles['window_size_x'], data_particles['window_size_y']]) / data_particles['window_size_x']
        # figsize=None

    fig.set_size_inches(*figsize)
        # ax.set_xlim(particles[:, 0].min(), particles[:, 0].max())
        # ax.set_ylim(particles[:, 1].min(), particles[:, 1].max())
    # if 'unwrap' in file:
    #     ax.set_xlim(-data_particles['window_size_x'], 2*data_particles['window_size_x'])
    #     ax.set_ylim(-data_particles['window_size_y'], 2*data_particles['window_size_y'])
    # else:
    #     ax.set_xlim(0, data_particles['window_size_x'])
    #     ax.set_ylim(0, data_particles['window_size_y'])
    # any ylim here is overridden later
        

        # print(frame.shape)
        # frame = frame.transpose([1, 0])
        # print(frame.shape)
    if not no_stack:
        frame = stack[0, :, :][::-1, :]
        print('^^^^^^ ahhh why is this here')

        preprocessing.stack_movie.show_single_frame(ax=ax, frame=frame, pixel_size=pixel_size,
                                                        window_size_x=data_particles['window_size_x'], window_size_y=data_particles['window_size_y'],
                                                        hide_scale_bar=True)

    print('adding outlines')
    add_outlines(0, ax)

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):

        fig, ax = plt.subplots(1, 1)

        go(file, ax, fig)

        if CROP:
            ax.set_xlim(X_START*pixel_size, (X_START+CROP)*pixel_size)
            ax.set_ylim(Y_START*pixel_size, (Y_START+CROP)*pixel_size)

        common.save_fig(fig, f'particle_detection/figures_png/frame1_{file}.png', dpi=300, only_plot=True)