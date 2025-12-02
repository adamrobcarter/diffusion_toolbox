import common
import numpy as np
import tqdm
import scipy.ndimage

if __name__ == '__main__':
    num_pixels = 2048 // 8
    pixel_size = 0.325
    L = num_pixels * pixel_size
    phi = 0.02
    sigma = 4
    radius = sigma/2
    time_step = 1
    D = common.stokes_einstein_D(sigma)
    num_timesteps = 500
    signal_to_noise = 0.01

    # for v in [0, 0.1, 0.2]:
    # for v in [0, 0.1, 0.2]:
    # for v in [0.04]:
    #     for signal_to_noise in [0.02]:
        # for signal_to_noise in [1]:
    for amplitude in [4]: # seems to be px
        for freq in [3.5]:

            num_particles = int(L**2 * 4 / np.pi * phi / sigma**2)

            rng = np.random.default_rng()

            print('getting data')
            stepsize = np.sqrt( 2 * D * time_step )
            print(f'stepsize {stepsize:.3g}')
            steps_x = rng.normal(0, stepsize, size=(num_particles, num_timesteps))
            startpoints_x = rng.uniform(0, L, size=(num_particles),              )
            steps_y = rng.normal(0, stepsize, size=(num_particles, num_timesteps))
            startpoints_y = rng.uniform(0, L, size=(num_particles),              )

            print('summing')
            particle_x = startpoints_x[:, np.newaxis] + np.cumsum(steps_x, axis=1)
            particle_y = startpoints_y[:, np.newaxis] + np.cumsum(steps_y, axis=1)

            print('modding into box')
            particle_x = particle_x % L
            particle_y = particle_y % L

            pixel_x = np.arange(0, num_pixels) * pixel_size
            pixel_y = np.arange(0, num_pixels) * pixel_size
            pixel_xx, pixel_yy = np.meshgrid(pixel_x, pixel_y)

            def single_gaussian_particle(x, y):
                # i = np.zeros([num_pixels, num_pixels])
                r2 = (x - pixel_xx)**2 + (y - pixel_yy)**2
                # i[r2 < 3*sigma] = np.exp(- r2[r2 < 3*sigma] / sigma)
                # return i
                i = np.exp(- r2 / radius)
                assert np.isfinite(i).all()
                return i
            
            def single_spherical_particle(x, y):
                r2 = (x - pixel_xx)**2 + (y - pixel_yy)**2
                r2[r2 > radius**2] = radius**2 # hack to prevent sqrt of negative number

                i = np.sqrt(radius**2 - r2)
                i /= radius # normalise so max is 1
                return i

            background = rng.uniform(0, 1/signal_to_noise, size=(int((num_pixels)*10), int((num_pixels+2*amplitude*num_timesteps)*10)))
            background = scipy.ndimage.uniform_filter(background, size=30, mode='reflect')
            assert np.isfinite(background).all()
            # ^ this is the true background which we reproject each time

            def background_at_t(t):
                # apply the background velocity
                offset = amplitude * ( np.sin(2*np.pi*freq*t) + 1 ) # + 1 cause the offset should be positive

                bkg = background[0:num_pixels*10, int(offset):int(num_pixels*10+offset)]
                # now we gotta do the mean per chunk
                # you could do this better with reshape or summat but i cba
                i = np.zeros([num_pixels, num_pixels])
                for x in range(num_pixels):
                    for y in range(num_pixels):
                        i[x, y] = bkg[x*10 : (x+1)*10, y*10 : (y+1)*10].mean()

                assert np.isfinite(i).all()
                return i

            stack = np.zeros([num_timesteps, num_pixels, num_pixels])
            progress = tqdm.tqdm(total=num_particles * num_timesteps, desc='adding background and particles')
            for t in range(num_timesteps):
                stack[t, :, :] += background_at_t(t)
                for p in range(num_particles):
                    stack[t, :, :] += single_spherical_particle(particle_x[p, t], particle_y[p, t])
                    progress.update()
            progress.close()

            # print(stack.mean(), stack.std(), stack.min(), stack.max()0)
            common.term_hist(stack)

            common.save_data(f'preprocessing/data/stack_sim_psiche_vib{freq}_a{amplitude}.npz',
                            stack=stack, pixel_size=pixel_size, time_step=time_step,
                            signal_to_noise=signal_to_noise,
                            particle_diameter=sigma)
            common.save_data(f'preprocessing/data/stack_sim_psiche_vib{freq}_a{amplitude}_small.npz',
                            stack=stack[::int(num_timesteps/50)], pixel_size=pixel_size, time_step=time_step,
                            signal_to_noise=signal_to_noise,
                            particle_diameter=sigma)