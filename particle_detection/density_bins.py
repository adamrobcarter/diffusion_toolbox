import common
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        time_column = data['dimension']
        t_max = particles[:, time_column].max()

        fig, ax = plt.subplots()

        bins = np.linspace(0, data['window_size_x'], num=50)

        for target_t in [0, t_max//3, t_max-1]:
            particles_at_t = particles[particles[:, time_column] == target_t, :]
            # counts, _ = np.histogram(particles_at_t[:, 0])

            ax.hist(particles_at_t[:, 0], bins=bins, alpha=0.5, label=f't={target_t}')

        ax.legend()
        common.save_fig(fig, f'particle_detection/figures_png/density_bins_{file}.png')