import common
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys
import warnings

def show_frame(fig, ax, stack, pixel_size, particles, radius, timestep, file):

    # crop = 0.18
    # crop = 1.0
    # particles_inside = ( particles[:, 0] < crop*stack.shape[1]*pixel_size ) & ( particles[:, 1] < crop*stack.shape[2]*pixel_size )
    # particles = particles[particles_inside, :]
    # radius    = radius   [particles_inside]
    # stack = stack[:, :int(stack.shape[1]*crop), :int(stack.shape[2]*crop)]

    # im = ax.imshow(stack[timestep, :, :], vmin=0, vmax=stack.max()/100)#, cmap=matplotlib.cm.coolwarm)#, vmin=-2, vmax=2)
    # there was a minus sign below for nice presentation
    # print(stack[timestep, :, :].min(), stack[timestep, :, :].max(), stack[timestep, :, :].mean())
    # print('mean', stack[timestep, :, :].mean())
    # assert 0.4 < stack[timestep, :, :].mean() < 0.5
    vmin = 0
    vmax = 1
    if file.startswith('marine'):
        pass
    elif file == 'pierre_exp':
        vmin = 0.3
        vmax = 0.6
    im = ax.imshow(stack[timestep, :, :], cmap=matplotlib.cm.Greys, vmin=vmin, vmax=vmax)#, vmin=-2, vmax=2)
    # im = ax.imshow(stack[timestep, :, :], cmap=matplotlib.cm.Greys, vmin=0, vmax=1)#, vmin=-2, vmax=2)
    # im = ax.imshow(stack[timestep, :, :], cmap=matplotlib.cm.Greys)#, vmin=-2, vmax=2)
    # fig.colorbar(im)

    TIME_INDEX = 2
    X_INDEX = 1
    Y_INDEX = 0 # why not 0/1?

    particles_at_t = particles[:, TIME_INDEX] == timestep
    # print(particles_at_t.sum())
    if particles_at_t.sum() == 0:
        warnings.warn('No particles found at timestep')
    # assert particles_at_t.sum() > 0
    # plt.scatter(particles[particles_at_t, X_INDEX]/pixel_size, particles[particles_at_t, Y_INDEX]/pixel_size, s=50*radius[particles_at_t]**2*pixel_size**2,
    #             facecolors='none', edgecolors='red', alpha=0.5, linewidth=0.8)
    
    x = particles[particles_at_t, X_INDEX]/pixel_size
    y = particles[particles_at_t, Y_INDEX]/pixel_size
    r = radius[particles_at_t] * np.sqrt(2)

    outline = True
    if outline:
        # print(r[0])
        # print(r[0, 0])
        # # print(r[0, 4])
        # print(r[4, 0])
        # print(r[4, 0]*1.5)
        circles = [plt.Circle((x[i],y[i]), radius=r[i]*2) for i in range(particles_at_t.sum())]
        # c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)
        c = matplotlib.collections.PatchCollection(circles, edgecolor='red', facecolor='none', linewidth=0.5, alpha=0.5)
    else:
        circles = [plt.Circle((x[i],y[i]), radius=r[i]) for i in range(particles_at_t.sum())]
        c = matplotlib.collections.PatchCollection(circles, facecolor='red', alpha=0.5)

    ax.add_collection(c)

    common.add_scale_bar(ax, pixel_size)

if __name__ == '__main__':
    for datasource in sys.argv[1:]:
        datasource2 = datasource
        if datasource.endswith('_trackpy'):
            datasource2 = datasource.split('_trackpy')[0]
        data = common.load(f'preprocessing/data/stack_{datasource2}.npz')
        stack             = data['stack']
        pixel_size        = data['pixel_size']
        particle_diameter = data['particle_diameter']
        num_timesteps = stack.shape[0]
        
        stack = stack - stack.mean(axis=0)
        stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

        # frame_average = stack.mean(axis=0)
        # stack = stack - frame_average[np.newaxis, :, :]
        # stack += frame_average / 10
        # stack = frame_average / 10 - stack
        
        # stack = stack - stack.mean(axis=0) # remove space background
        
        # stack = np.interp(stack, (stack.min(), stack.max()), (0, 1)) # convert to 0->1 range

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

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        # fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        fig.tight_layout()

        # for timestep in range(0, 1000, 10):
        if True:
            timestep = num_timesteps // 2

            show_frame(fig, ax, stack, pixel_size, particles, radius, timestep)

            # plt.axis('off') # hides all axes, leaving just figure
            fig.savefig(f'particle_detection/figures_png/particles_{datasource}.png', dpi=300, bbox_inches='tight', pad_inches=0)