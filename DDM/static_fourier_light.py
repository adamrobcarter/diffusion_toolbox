import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize

FIRST_FRAME = False
FRAME_DIFF = True
REMOVE_BKG = False

def do_static_fourier():

    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data['stack']
        pixel_size = data['pixel_size']

    if REMOVE_BKG:
        print('subtracting mean')
        stack = stack - stack.mean(axis=0)

    if FIRST_FRAME:
        if FRAME_DIFF:
            images = stack[[1], :, :] - stack[[0], :, :]
        else:
            images = stack[[0], :, :]
    else:
        if FRAME_DIFF:
            images = stack[1::5, :, :] - stack[:-1:5, :, :]
        else:
            images = stack[:, :, :]
    del stack
    num_frames = images.shape[0]
    
    # image = image[::4, ::4]
    # pixel_size *= 4

    # rng = np.random.default_rng()
    # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
    #                     2 * np.pi * np.arange(200) / 34)
    # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

    print(images.shape)
    fx, fy, fouriers = common.fourier_2D(images, pixel_size, (1, 2))
    del images
    # print(fouriers.shape)
    fourier = fouriers.mean(axis=0)
    # fx = fx.mean(axis=0)
    # fy = fy.mean(axis=0)
    print('f', fourier.shape, fx.shape)


    title = common.name(file)
    title += f'\n{num_frames} frames'
    title += ', diff' if FRAME_DIFF else ', nodiff'
    title += ', bkg rem' if REMOVE_BKG else ', no bkg rem'

    # a = np.abs(fourier)
    import scipy.fft
    a = np.abs(scipy.fft.fftshift(fourier))**2
    # TODO: should we be doing fftshift in the DDM code?
    # print(a.mean(), a.std())
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    # ax.semilogv()
    log_a = np.log(a)
    # im = ax.imshow(log_a, vmin=log_a.mean()-2*log_a.std(), vmax=log_a.mean()+2*log_a.std(), extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
    im = ax.imshow(log_a, extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
    fig.colorbar(im)
    # plt.imshow(image)
    ax.set_title(title)
    # ax.set_xlabel('$k_x$')
    # ax.set_ylabel('$k_y$')
    common.save_fig(fig, f'DDM/figures_png/static_fourier_{file}.png')

    
    # radial average
    f = scipy.fft.fftshift( np.sqrt( fx**2 + fy**2 ) )
    print('f', f.min(), f.max())
    bins = np.linspace(0, f.max()/np.sqrt(2), 1000)[1:]
    print('bin_sep', 2*np.pi/(bins[1]-bins[0]))
    ff = np.abs(fourier)**2
    print('f, a', f.shape, a.shape)
    v, bins, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), bins=bins)
    err, _, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), statistic='std', bins=bins)
    n, _, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), statistic='count', bins=bins)
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    x = (bins[1:] + bins[:-1]) / 2
    k = 2 * np.pi * x
    ax.scatter(k, v, s=2)
    ax.errorbar(k, v, yerr=err/np.sqrt(n), linestyle='none')
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel(r'$k$ ($\mathrm{\mu m}^{-1}$)')
    ax.set_ylabel(r'$\langle I(k)^* I(k) \rangle$')

    # ax.vlines(2*np.pi/1.5, *ax.get_ylim())

    if particle_diameter:
        print(particle_diameter)
        print('onit', 2*np.pi/particle_diameter)
        ax.vlines(2*np.pi/particle_diameter, v.min(), v.max())
    # ax.vlines(2*np.pi/data['pixel_size'], v.min(), v.max())
    ax.set_ylim(v.min()/1.1, v.max()*1.1)

    ax.grid()

    end = 2e-1
    end_index = np.argmax(k > end)
    print(end_index)
    # func = lambda x, a, b, c: x**a
    # popt, pcov = scipy.optimize.curve_fit(func, k[:end_index], v[:end_index])
    # plt.plot(k[:end_index], func(k[:end_index], *popt), color='black', label=f'k**{popt[0]:.3f}')
    
    # func = lambda x, b, a: b * x + a
    # popt, pcov = scipy.optimize.curve_fit(func, np.log10(k[:end_index]), np.log10(v[:end_index]))
    # plt.plot(k[:end_index], 10**func(np.log10(k[:end_index]), *popt), color='black', label=f'k**{popt[0]:.3f}')
    
    
    plt.legend()

    realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
    # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
    realspace_ax.set_xlabel(r'$2\pi/k$ ($\mathrm{\mu m}$)')

    ax.set_xlim(k.min()/1.1, max(k.max(), 2*np.pi/pixel_size)*1.1)
    
    ax.set_title(title)
    # average_I_sq = image.mean()**2
    # print(average_I_sq, '<I^2>')

    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/static_fourier_av_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'DDM/figures_png/static_fourier_av_{file}.png', dpi=200)

if __name__ == '__main__':

    do_static_fourier()