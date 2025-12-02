import common
import numpy as np
import matplotlib.pyplot as plt

NUM_FRAMES = int(1e10)

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1)

    for file_i, file in enumerate(files := common.files_from_argv('particle_linking/data', 'trajs_')):
        if '_div' in file:
            continue
        if 'theta0.66' in file:
            continue

        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particles = data['particles']
        theta     = data['theta']

        TIME_COLUMN = data.get('dimension', 2)
        ID_COLUMN   = data.get('dimension', 2) + 1

        assert len(np.unique(particles[:, ID_COLUMN])) == 1, 'more than one particle found'

        ax.plot(particles[:NUM_FRAMES, TIME_COLUMN], 2*particles[:NUM_FRAMES, 2]/data['particle_diameter'], label=fr'$\theta={theta}$',
                color=common.colormap(file_i, 0, len(files)))

    ax.set_ylabel('$z / a$')
    ax.set_xlabel('$t$')
    ax.legend(fontsize=8)
    common.save_fig(fig, f'particle_linking/figures_png/z_{file}.png', dpi=300)
