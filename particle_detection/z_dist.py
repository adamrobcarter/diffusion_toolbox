import common
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    a = data['particle_diameter'] / 2
    z = particles[:, 2]
    zoa = z/a

    label = f'{file}\n$z = {z.mean():.3f} \\pm {z.std():.3f} = ({zoa.mean():.3f} \\pm {zoa.std():.3f})a$'
    ax.hist(zoa, bins=50, label=label, alpha=0.7, density=True)

ax.set_xlabel('$z/a$')
ax.legend()
common.save_fig(fig, f'particle_detection/figures_png/z_dist_{file}.png')