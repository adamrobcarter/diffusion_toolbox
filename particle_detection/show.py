import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys

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

    density = particles.shape[0]/num_timesteps / (stack.shape[0]*stack.shape[1]*pixel_size**2)
    pack_frac = np.pi/4 * density * particle_diameter**2
    print(f'packing fraction phi = {pack_frac:.3f}')

    plt.figure(figsize=(12, 12))
    plt.tight_layout()

    # for timestep in range(0, 1000, 10):
    if True:
        timestep = num_timesteps // 2

        im = plt.imshow(stack[timestep, :, :], vmin=0, vmax=stack.max()/8)#, cmap=matplotlib.cm.coolwarm)#, vmin=-2, vmax=2)
        # plt.savefig('figures_png/start.png')
        plt.colorbar(im)

        TIME_INDEX = 2
        X_INDEX = 1
        Y_INDEX = 0 # why not 0/1?

        particles_at_t = particles[:, TIME_INDEX] == timestep
        plt.scatter(particles[particles_at_t, X_INDEX]/pixel_size, particles[particles_at_t, Y_INDEX]/pixel_size, s=50*radius[particles_at_t]**2*pixel_size**2,
                    facecolors='none', edgecolors='red', alpha=0.5, linewidth=0.5)
        
        plt.savefig(f'particle_detection/figures_png/particles_{datasource}.png', dpi=300)