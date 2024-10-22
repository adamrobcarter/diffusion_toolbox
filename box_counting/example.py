import common
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.patches

def go(file, box, sep, export_destination=None, label_boxes=True):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']

    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data['stack']
    pixel_size = data['pixel_size']

    FRAME = 0

    fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    ax.imshow(stack[0, :, :], cmap=matplotlib.cm.Greys, interpolation='none')
    # plt.imshow(stack.min(axis=0))
    
    # excess = stack - stack.min(axis=0)
    # print(stack[:, :, timestep].mean())


    # common.add_scale_bar(ax, data['pixel_size'], color='black')

    box_size = box / pixel_size
    sep_size = sep / pixel_size

    window_size = 2 * 36 + 2.5 * 10

    

    # box_size = 30 / pixel_size
    # box_size = 22 / pixel_size
    # # sep_size = 20 / pixel_size
    # sep_size = 2.88 * 4 / pixel_size

    window_size = 356

    particles_at_t = particles[particles[:, 2] == FRAME]

    x = sep_size*0.8
    while x < window_size:
        y = sep_size*0.8
        while y < window_size:

            N = 0
            for particle in particles_at_t:
                if x < particle[1]/pixel_size < x+box_size and y < particle[0]/pixel_size < y+box_size:
                    N += 1

            color = matplotlib.cm.afmhot(0.3)
            # color = 'xkcd:dark orange'
            rect = matplotlib.patches.Rectangle((x, y), box_size, box_size, edgecolor=color, facecolor='none')
            ax.add_patch(rect)

            if label_boxes:
                ax.text(x+box_size/2, y+box_size/2, f'$N(t) = {N}$', fontsize=9, color=color, horizontalalignment='center', verticalalignment='center')
                ax.text(x+box_size/2, y+box_size/8, fr'$L = {box_size*pixel_size:.0f}\mathrm{{\mu m}}$', fontsize=7, color=color, horizontalalignment='center', verticalalignment='baseline')

            y += box_size + sep_size
        x += box_size + sep_size

    ax.set_ylim(000, window_size)
    ax.set_xlim(000, window_size)

    if export_destination:
        common.save_fig(fig, export_destination, only_plot=True, hide_metadata=True)
    common.save_fig(fig, f'box_counting/figures_png/example_{file}.png', dpi=300, only_plot=True)
    
if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file, 76, 40)