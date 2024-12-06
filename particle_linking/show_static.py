import common
# import preprocessing.stack_movie
# import particle_detection.show_movie
import sys
import matplotlib.pyplot as plt
import tqdm
import matplotlib.cm

MAX_TIME = 1000
CROP = 150
X_START = 250
Y_START = 50

def go(file, destination_trajs, destination_positions, destination_intensities):
    data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data_stack['stack']
    pixel_size = data_stack['pixel_size']
    #     found_stack = True
    # except FileNotFoundError:
    #     found_stack = False

    data_particles = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data_particles['particles']

    # crop
    # stack = stack[:, X_START:X_START+CROP, Y_START:Y_START+CROP]

    stack = stack - stack.mean(axis=0) # remove space background
    
    # def add_outlines(timestep, ax):
    #     particle_detection.show.add_particle_outlines(ax, pixel_size, particles, radius, timestep)

    filename = f'movie_linked_{file}'

    fig, ax = plt.subplots(1, 1)

    ax.imshow(stack[0, :, :], cmap=matplotlib.cm.Greys)
    if destination_intensities:
        common.save_fig(fig, destination_intensities, only_plot=True)

    circles = []

    # find particles in view
    ids = []
    for row_i in tqdm.trange(particles.shape[0]):
        row = particles[row_i, :]
        if row[2] == 0:
            # if X_START*pixel_size < row[0] < (X_START+CROP)*pixel_size:
            #     if Y_START*pixel_size < row[1] < (Y_START+CROP)*pixel_size:
                    ids.append(row[3])

                    s = ax.scatter(row[1]/pixel_size, row[0]/pixel_size, edgecolors='red', facecolors='none', s=1500)

                    circles.append(s)

    if destination_positions:
        common.save_fig(fig, destination_positions, only_plot=True)

    for s in circles:
        s.remove()

    # find trajs
    for id in tqdm.tqdm(ids):
        particles_this_id = particles[particles[:, 3] == id, :]
        particles_this_time = particles_this_id[particles_this_id[:, 2] < MAX_TIME, :]
        x = particles_this_time[:, 0] / pixel_size - X_START
        y = particles_this_time[:, 1] / pixel_size - Y_START

        ax.plot(y, x)
        ax.plot(x, y)

    ax.set_xlim(X_START, X_START+CROP)
    ax.set_ylim(Y_START, Y_START+CROP)

    common.save_fig(fig, destination_trajs, only_plot=True)

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        destination = f'particle_linking/figures_png/linked_static_{file}.png'
        go(file,
           f'particle_linking/figures_png/static_linked_{file}.png',
           f'particle_linking/figures_png/static_positi_{file}.png',
           f'particle_linking/figures_png/static_intens_{file}.png',
        )