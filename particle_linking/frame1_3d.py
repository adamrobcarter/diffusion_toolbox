import common
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

POINTS_ON_SPHERE = 10

def plot_sphere(ax, coord, r):
    u = np.linspace(0, 2 * np.pi, POINTS_ON_SPHERE)
    v = np.linspace(0, np.pi, POINTS_ON_SPHERE//2)

    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones_like(u), np.cos(v))

    x += coord[0]
    y += coord[1]
    z += coord[2]

    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')

for file in common.files_from_argv('particle_linking/data', 'trajs_'):
    data = common.load(f'particle_linking/data/trajs_{file}.npz')

    assert data.get('dimension') == 3
    time_column = 3
    id_columns = 4

    particles = data['particles']

    fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
    ax = fig.add_subplot(projection='3d')

    particles_t0 = particles[particles[:, time_column] == 0, :]

    for row_i in range(particles_t0.shape[0]):
        plot_sphere(ax, particles_t0[row_i, [0, 1, 2]], r=data['particle_diameter']/2)

    common.save_fig(fig, f'particle_linking/figures_png/frame1_{file}.png')