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

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')

        assert data.get('dimension') == 3
        time_column = 3
        id_columns = 4

        particles = data['particles']
        num_timesteps = len(np.unique(particles[:, time_column]))
        timestep = np.unique(particles[:, time_column])[1]
        print('timestep', timestep, 'num_timesteps', num_timesteps)
        particles[:, time_column] /= timestep
        particles[:, time_column] = np.round(particles[:, time_column]) # b/c of irrational timesteps

        fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
        ax = fig.add_subplot(projection='3d')
        ax.dist = 5 # zooms in. 10 is the default. smaller == more zoom

        def show_frame(frame):
            ax.clear()
            particles_t = particles[particles[:, time_column] == frame, :]

            for row_i in range(particles_t.shape[0]):
                plot_sphere(ax, particles_t[row_i, [0, 1, 2]], r=data['particle_diameter']/2)

            ax.text2D(0.05, 0.95, f't={timestep*frame:.0f}s', transform=ax.transAxes)

        common.save_gif(show_frame, range(0, num_timesteps, 100), fig, f'particle_linking/figures_png/movie_{file}.gif', fps=10)