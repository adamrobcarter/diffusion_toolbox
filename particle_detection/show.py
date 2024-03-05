import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys

def show_frame(fig, ax, stack, pixel_size, particles, radius, timestep):
    im = ax.imshow(stack[timestep, :, :], vmin=0, vmax=stack.max()/2)#, cmap=matplotlib.cm.coolwarm)#, vmin=-2, vmax=2)
    # fig.colorbar(im)

    TIME_INDEX = 2
    X_INDEX = 1
    Y_INDEX = 0 # why not 0/1?

    particles_at_t = particles[:, TIME_INDEX] == timestep
    # plt.scatter(particles[particles_at_t, X_INDEX]/pixel_size, particles[particles_at_t, Y_INDEX]/pixel_size, s=50*radius[particles_at_t]**2*pixel_size**2,
    #             facecolors='none', edgecolors='red', alpha=0.5, linewidth=0.5)
    x = particles[particles_at_t, X_INDEX]/pixel_size
    y = particles[particles_at_t, Y_INDEX]/pixel_size
    r = radius[particles_at_t] * np.sqrt(2)
    circles = [plt.Circle((x[i],y[i]), radius=r[i]) for i in range(particles_at_t.sum())]
    c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)
    ax.add_collection(c)

if __name__ == '__main__':
    for datasource in sys.argv[1:]:
        data = common.load(f'preprocessing/data/stack_{datasource}.npz')
        stack             = data['stack']
        pixel_size        = data['pixel_size']
        particle_diameter = data['particle_diameter']
        num_timesteps = stack.shape[0]

        frame_average = stack.mean(axis=0)
        stack = stack - frame_average[np.newaxis, :, :]

        print('stack max', stack.max())

        data = common.load(f'particle_detection/data/particles_{datasource}.npz')
        particles = data['particles']
        radius = data['radius']
        print(particles.shape)

        ##### radii are wrong! pls fix!!

        density = particles.shape[0]/num_timesteps / (stack.shape[0]*stack.shape[1]*pixel_size**2)
        pack_frac = np.pi/4 * density * particle_diameter**2
        print(f'packing fraction phi = {pack_frac:.3f}')

        print(f'radius min', radius.min(), 'median', np.median(radius), 'max', radius.max())
        print(f'radius min', radius.min()* np.sqrt(2), 'median', np.median(radius)* np.sqrt(2), 'max', radius.max()* np.sqrt(2))

        # fig, ax = plt.subplots(1, 1, figsize=(12, 12))
        fig, ax = plt.subplots(1, 1)
        fig.tight_layout()

        # for timestep in range(0, 1000, 10):
        if True:
            timestep = num_timesteps // 2

            show_frame(fig, ax, stack, pixel_size, particles, radius, timestep)

            fig.savefig(f'particle_detection/figures_png/particles_{datasource}.png', dpi=1200)