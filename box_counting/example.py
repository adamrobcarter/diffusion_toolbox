import common
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.patches

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']

    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    pixel_size = data['pixel_size']

    FRAME = 0

    fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    ax.imshow(stack[0, :, :], cmap=matplotlib.cm.Greys)
    # plt.imshow(stack.min(axis=0))
    
    # excess = stack - stack.min(axis=0)
    # print(stack[:, :, timestep].mean())

    window_size = 600
    ax.set_ylim(000, window_size)
    ax.set_xlim(000, window_size)

    # common.add_scale_bar(ax, data['pixel_size'], color='black')


    box_size = 60 / pixel_size
    sep_size = 20 / pixel_size

    particles_at_t = particles[particles[:, 2] == FRAME]

    x = sep_size*0.8
    while x < window_size:
        y = sep_size*0.8
        while y < window_size:

            N = 0
            for particle in particles_at_t:
                if x < particle[1]/pixel_size < x+box_size and y < particle[0]/pixel_size < y+box_size:
                    N += 1

            rect = matplotlib.patches.Rectangle((x, y), box_size, box_size, edgecolor='blue', facecolor='none', alpha=0.7)
            ax.add_patch(rect)
            ax.text(x+box_size/2, y+box_size/2, f'$N(t) = {N}$', fontsize=8, color='blue', alpha=0.7, horizontalalignment='center', verticalalignment='center')
            ax.text(x+box_size/2, y+10, fr'$L = {box_size*pixel_size:.0f}\mathrm{{\mu m}}$', fontsize=6, color='blue', alpha=0.7, horizontalalignment='center', verticalalignment='baseline')

            y += box_size + sep_size
        x += box_size + sep_size


    common.save_fig(fig, f'box_counting/figures_png/example_{file}.png', dpi=300, only_plot=True)
    common.save_fig(fig, f'/home/acarter/presentations/intcha24/figures/boxcounting_example_{file}.png', dpi=300, only_plot=True, hide_metadata=True)
    