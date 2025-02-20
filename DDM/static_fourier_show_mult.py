import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

SAVE_TWO_D_PLOT = True

def form_factor_sphere(q, R):
    Fs = 3 * (np.sin(q*R) - q*R*np.cos(q*R)) / (q*R)**3
    return Fs**2

def show_static_fourier(file, ax, color):

    data = common.load(f'DDM/data/static_fourier_{file}.npz')
    pixel_size        = data['pixel_size']
    num_frames        = data['num_frames_used']
    particle_diameter = data.get('particle_diameter')
    fourier           = data['fourier']
    k_x                = data['fx'] * 2 * np.pi
    k_y                = data['fy'] * 2 * np.pi

    assert np.isfinite(fourier).all()

    fourier = np.fft.fftshift(fourier)
    k_x = np.fft.fftshift(k_x)
    k_y = np.fft.fftshift(k_y)

    title = common.name(file)
    title += f', {num_frames} frames'
    # title += ', diff' if FRAME_DIFF else ', nodiff'
    # title += ', bkg rem' if REMOVE_BKG else ', no bkg rem'

    # a = np.abs(fourier)
    fourier_abs = np.abs(fourier)**2
    # TODO: should we be doing fftshift in the DDM code?
    # print(a.mean(), a.std())
    # fig, ax = plt.subplots(1, 1, figsize=(5, 5)) # this is for the 2D plot
    # ax.semilogv()
    log_fourier_abs = np.log(fourier_abs)
    # im = ax.imshow(log_a, vmin=log_a.mean()-2*log_a.std(), vmax=log_a.mean()+2*log_a.std(), extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
    # if SAVE_TWO_D_PLOT:
    #     # im = ax.imshow(log_a, extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
    #     im = ax.pcolormesh(k_x, k_y, log_fourier_abs, shading='nearest')
    #     fig.colorbar(im)
    #     ax.set_title(title)
    #     # ax.set_xlabel('$k_x$')
    #     # ax.set_ylabel('$k_y$')
    #     filename = f'static_fourier_{file}'
    #     common.save_fig(fig, f'DDM/figures_png/{filename}.png')

    ##################################
    # radial average
    ##################################
    title = common.name(file)
    # title += f'\n{num_frames} frames'
    title += f'\n'
    
    # f = scipy.fft.fftshift( np.sqrt( k_x**2 + fy**2 ) )
    k = np.sqrt(k_x**2 + k_y**2)

    min_k = np.abs(k_x[0, 1] - k_x[0, 0])

    print('f', k.min(), k.max())
    bins = np.arange(min_k, k.max()/np.sqrt(2), step=min_k)
    print('bin_sep', 2*np.pi/(bins[1]-bins[0]))
    # ff = np.abs(fourier)**2
    print('f, a', k.shape, fourier_abs.shape)
    k_flat = k.flatten()
    fourier_abs_flat = fourier_abs.flatten()
    
    v, bins, _ = scipy.stats.binned_statistic(k_flat, fourier_abs_flat, bins=bins)
    assert len(bins.shape) == 1
    print(v)
    assert np.isfinite(v).all()
    err, _, _ = scipy.stats.binned_statistic(k.flatten(), fourier_abs.flatten(), statistic='std', bins=bins)
    n, _, _ = scipy.stats.binned_statistic(k.flatten(), fourier_abs.flatten(), statistic='count', bins=bins)
    k_bin_mids = (bins[1:] + bins[:-1]) / 2
    # k = 2 * np.pi * k_bin_mids
    ax.scatter(k_bin_mids, v, s=2, color=color, label=file)
    ax.errorbar(k_bin_mids, v, yerr=err/np.sqrt(n), linestyle='none', color=color)
    
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel(r'$k$ ($\mathrm{\mu m}^{-1}$)')
    ax.set_ylabel(r'$\langle I(k)^* I(k) \rangle$')

    # if particle_diameter:
    #     print(particle_diameter)
    #     print('onit', 2*np.pi/particle_diameter)
    #     ax.vlines(2*np.pi/particle_diameter, v.min(), v.max())
    # ax.vlines(2*np.pi/data['pixel_size'], v.min(), v.max())
    ax.set_ylim(v.min()/1.2, v.max()*1.2)

    ax.grid(alpha=0.3)

    end = 2e-1
    end_index = np.argmax(k > end)
    print(end_index)
    # func = lambda x, a, b, c: x**a
    # popt, pcov = scipy.optimize.curve_fit(func, k[:end_index], v[:end_index])
    # plt.plot(k[:end_index], func(k[:end_index], *popt), color='black', label=f'k**{popt[0]:.3f}')
    
    # func = lambda x, b, a: b * x + a
    # popt, pcov = scipy.optimize.curve_fit(func, np.log10(k[:end_index]), np.log10(v[:end_index]))
    # plt.plot(k[:end_index], 10**func(np.log10(k[:end_index]), *popt), color='black', label=f'k**{popt[0]:.3f}')
    

    realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
    # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
    realspace_ax.set_xlabel(r'$2\pi/k$ ($\mathrm{\mu m}$)')

    ax.set_xlim(k.min()/1.1, max(k.max(), 2*np.pi/pixel_size)*1.1)
    
    ax.set_title(title + str(data.get('NAME', '')))
    # average_I_sq = image.mean()**2
    # print(average_I_sq, '<I^2>')

    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/static_fourier_av_{file}.pdf', hide_metadata=True)
    return particle_diameter, k_bin_mids, v # dirty hack

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    
    for i, file in enumerate(files := common.files_from_argv('DDM/data', 'static_fourier_')):
        color = common.colormap(i, 0, len(files))

        particle_diameter, k_bin_mids, v = show_static_fourier(file, ax, color)

    
        
    if particle_diameter:
        Pq = form_factor_sphere(k_bin_mids, particle_diameter/2)
        index = np.argmax(k_bin_mids > 1e0)
        Pq *= v[index] / Pq[index]
        ax.plot(k_bin_mids, Pq, label=f'$P(q)$ ($2R={particle_diameter:.1f}\mathrm{{\mu m}}$)')
        print(form_factor_sphere(k_bin_mids, 4))

    common.add_exponential_index_indicator(ax, exponent=-4, anchor=(1, 1), xlabel='k')
    
    plt.legend(fontsize=9)

    filename = f'static_fourier_av_{file}'
    common.save_fig(fig, f'DDM/figures_png/{filename}.png', dpi=200)