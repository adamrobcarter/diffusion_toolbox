import numpy as np
import common
import matplotlib.pyplot as plt

for file in common.files_from_argv('particle_detection/data/', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']

    hist, bins = np.histogram(particles[:, 2], bins=range(0, int(particles[:, 2].max())+2))

    fig, ax = plt.subplots()
    ax.plot(hist)
    ax.set_xlabel('time (frames)')
    ax.set_ylabel('number of particles')
    common.save_fig(fig, f'particle_detection/figures_png/totalnumber_{file}.png')