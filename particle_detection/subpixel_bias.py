import numpy as np
import matplotlib.pyplot as plt
import common

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles  = data['particles'] # rows of x,y,t
    pixel_size = data['pixel_size']
    particles[[0, 1], :] /= pixel_size # go back to pixel units not real space units
    fig, (ax_x, ax_y) = plt.subplots(1, 2, figsize=(6, 3))
    ax_x.hist(particles[:, 0] % 1)
    ax_y.hist(particles[:, 1] % 1)
    ax_x.set_title('x')
    ax_y.set_title('y')
    ax_x.set_xlim(0, 1)
    ax_y.set_xlim(0, 1)
    fig.suptitle(f'subpixel bias {file}')
    common.save_fig(fig, f'particle_detection/figures_png/subpixel_bias_{file}.png')