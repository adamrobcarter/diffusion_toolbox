import common
import numpy as np
import matplotlib.pyplot as plt

def go(file):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')
    particles = data['particles']

    assert data['dimension'] == 3
    time_column = 3
    id_column = 4
    assert np.all(particles[:, id_column] == 0)

    fig, ax = plt.subplots(1, 1)
    t = particles[:, time_column]
    for i in range(3):
        axis = 'xyz'[i]
        ax.plot(t, particles[:, 5+i], label=fr'$\phi_{axis}$')

    ax.legend()

    # ax.set_ylim(-10, 10)
    # ax.set_xlim(0, 1000)
        
    common.save_fig(fig, f'particle_linking/figures_png/rotation_{file}.png', dpi=300)
    

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        go(file)