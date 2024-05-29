import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']

    image = stack[0, :, :]
    if stack.shape[0] > 1:
        print('subbing back')
        image = image - stack.mean(axis=0)
    
    # image = image[::4, ::4]
    # pixel_size *= 4

    # rng = np.random.default_rng()
    # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
    #                     2 * np.pi * np.arange(200) / 34)
    # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

    print(image.shape)
    fx, fy, fourier = common.fourier_2D(image, pixel_size, (0, 1))



    # a = np.abs(fourier)
    import scipy.fft
    a = np.abs(scipy.fft.fftshift(fourier))**2
    # TODO: should we be doing fftshift in the DDM code?
    # print(a.mean(), a.std())
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    # ax.semilogv()
    im = ax.imshow(np.log(a), extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
    plt.colorbar(im)
    # plt.imshow(image)
    common.save_fig(fig, f'DDM/figures_png/static_fourier_{file}.png')


    # radial average
    f = scipy.fft.fftshift( np.sqrt( fx**2 + fy**2 ) )
    bins = np.linspace(0, f.max()/np.sqrt(2), 1000)[1:]
    ff = np.abs(fourier)**2
    v, bins, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), bins=bins)
    err, _, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), statistic='std', bins=bins)
    n, _, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), statistic='count', bins=bins)
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    x = (bins[1:] + bins[:-1]) / 2
    k = 2 * np.pi * x
    ax.scatter(k, v)
    ax.errorbar(k, v, yerr=err/np.sqrt(n), linestyle='none')
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel(r'$k$ ($\mathrm{\mu m}^{-1}$)')
    ax.set_ylabel(r'$\langle I(k)^* I(k) \rangle$')



    if particle_diameter := data.get('particle_diameter'):
        print(particle_diameter)
        print('onit', 2*np.pi/particle_diameter)
        ax.vlines(2*np.pi/particle_diameter, v.min(), v.max())
    # ax.vlines(2*np.pi/data['pixel_size'], v.min(), v.max())
    ax.set_ylim(v.min()/1.1, v.max()*1.1)

    realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
    # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
    realspace_ax.set_xlabel(r'$2\pi/k$ ($\mathrm{\mu m}$)')

    ax.set_xlim(k.min()/1.1, max(k.max(), 2*np.pi/data['pixel_size'])*1.1)

    average_I_sq = image.mean()**2
    print(average_I_sq, '<I^2>')

    common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/static_fourier_av_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'DDM/figures_png/static_fourier_av_{file}.png')