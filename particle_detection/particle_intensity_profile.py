import common
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import tqdm

max_radius = 5
r_step_px = 1/2

for file in common.files_from_argv('particle_detection/data/', 'particles_'):
    data_particles = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data_particles['particles']
    
    data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data_stack['stack']
    x_pixels, y_pixels = np.meshgrid(range(stack.shape[2]), range(stack.shape[1]))
    #                                                  ^ again idk why this needs to be switched
    pixel_size = data_stack['pixel_size']
    
    assert stack[0, :, :].shape == x_pixels.shape, f'{stack[0, :, :].shape} != {x_pixels.shape}'
    assert stack[0, :, :].shape == y_pixels.shape, f'{stack[0, :, :].shape} != {y_pixels.shape}'

    max_radius_px = int(max_radius / pixel_size)
    print('max_radius_px', max_radius_px)
    # r_bins_px = np.array(range(0, max_radius_px, r_step_px))
    r_bins_px = np.arange(0, max_radius_px, r_step_px)

    max_k = 1 / (np.pi * max_radius)

    lines_to_use = range(0, particles.shape[0], 20)
    print(f'using {len(lines_to_use)} image sections')
    intensity_profiles = np.full((len(lines_to_use), len(r_bins_px)-1), np.nan)
    used_images = np.full((len(lines_to_use), max_radius_px*2, max_radius_px*2), np.nan)

    for index, line in tqdm.tqdm(enumerate(lines_to_use), total=len(lines_to_use)):
        line = particles[line, :]
        particle_x, particle_y, t = line[0], line[1], int(line[2])

        particle_x_px = particle_x / pixel_size
        particle_y_px = particle_y / pixel_size

        if (
            particle_x_px < max_radius_px or particle_x_px > stack.shape[1]-max_radius_px or
            particle_y_px < max_radius_px or particle_y_px > stack.shape[2]-max_radius_px
        ):
            continue # skip particles close to the edge for now

        image_to_use    = stack[t, int(particle_x_px-max_radius_px):int(particle_x_px+max_radius_px), int(particle_y_px-max_radius_px):int(particle_y_px+max_radius_px)]
        assert image_to_use.shape == (max_radius_px*2, max_radius_px*2), f'image_to_use.shape = {image_to_use.shape}. should be ({max_radius_px*2}, {max_radius_px*2})'
        x_pixels_to_use = y_pixels[int(particle_x_px-max_radius_px):int(particle_x_px+max_radius_px), int(particle_y_px-max_radius_px):int(particle_y_px+max_radius_px)]
        y_pixels_to_use = x_pixels[int(particle_x_px-max_radius_px):int(particle_x_px+max_radius_px), int(particle_y_px-max_radius_px):int(particle_y_px+max_radius_px)]
        #                 ^ i do not know why these are switched ffs ahaha
        assert x_pixels_to_use.shape == (max_radius_px*2, max_radius_px*2), f'x_pixels_to_use.shape = {x_pixels_to_use.shape}. should be ({max_radius_px*2}, {max_radius_px*2})'
        assert y_pixels_to_use.shape == (max_radius_px*2, max_radius_px*2), f'y_pixels_to_use.shape = {y_pixels_to_use.shape}. should be ({max_radius_px*2}, {max_radius_px*2})'

        # print(x_pixels_to_use.min(), particle_x_px, x_pixels_to_use.max())
        # print(y_pixels_to_use.min(), particle_y_px, y_pixels_to_use.max())
        # print()

        # find th distance between each pixel and the centre of the particle
        d_px = np.sqrt( (x_pixels_to_use-particle_x_px)**2 + (y_pixels_to_use-particle_y_px)**2)
        # print(d_px)

        # azimuthal average
        I_binned, d_binned, _ = scipy.stats.binned_statistic(d_px.flatten(), image_to_use.flatten(), bins=r_bins_px)
        # print(I_binned.shape, d_binned.sh)
        intensity_profiles[index, :] = I_binned

        used_images[index, :, :] = image_to_use

    fig, ((r_ax_2d, r_ax_1d), (k_ax_2d, k_ax_1d)) = plt.subplots(2, 2, figsize=(6, 6))
    avg_image = np.nanmean(used_images, axis=0)
    # note: don't just do binning on the average image, because you have lost subpixel accuracy if you do that
    r_ax_2d.imshow(avg_image)
    # ax_2d.imshow(d_px)
    r_bins_px_av = (r_bins_px[:-1] + r_bins_px[1:]) / 2 # middle of bin
    # print()
    intensity_profile = np.nanmean(intensity_profiles, axis=0)
    intensity_profile -= stack.mean()
    r_ax_1d.scatter(r_bins_px_av, intensity_profile)
    r_ax_1d.hlines(0, r_bins_px_av[0], r_bins_px_av[-1], zorder=-1, color='black', alpha=0.2)

    k_xx, k_yy, k_image = common.fourier_2D(avg_image, pixel_size, axes=(0, 1))
    print('k image.shape', k_image.shape)
    k_ax_2d.imshow(np.abs(k_image))

    test_f, test_fft = common.fourier(r_bins_px_av, intensity_profiles[0, :])
    k_profiles = np.full((len(lines_to_use), test_fft.size), np.nan)
    for i in range(len(lines_to_use)):
        print(r_bins_px_av.shape, intensity_profiles[i, :].shape, common.fourier(r_bins_px_av, intensity_profiles[i, :])[1].shape)
        k, k_profiles[i, :] = common.fourier(r_bins_px_av, intensity_profiles[i, :])
        # for sure there is a one function call that would work here instead of the loop

    k_profile = np.nanmean(k_profiles, axis=0)
    k_ax_1d.scatter(k, k_profile)
    k_ax_1d.semilogy()

    plt.savefig(f'particle_detection/figures_png/intensity_profile_{file}.png')

    np.savez(f'particle_detection/data/particle_intensity_profile_{file}.npz',
             profile=k_profile, k=k,
             max_radius=max_radius, r_step_px=r_step_px)
    



        