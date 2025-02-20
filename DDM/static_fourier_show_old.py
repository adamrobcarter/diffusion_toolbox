import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

SAVE_TWO_D_PLOT = False

def form_factor_sphere(q, R):
    Fs = 3 * (np.sin(q*R) - q*R*np.cos(q*R)) / (q*R)**3
    return Fs**2

def show_static_fourier(file):

    data = common.load(f'DDM/data/static_fourier_{file}.npz')
    pixel_size        = data['pixel_size']
    num_frames        = data['num_frames_used']
    particle_diameter = data.get('particle_diameter')
    fourier           = data['fourier']
    fx                = data['fx']
    fy                = data['fy']

    title = common.name(file)
    title += f', {num_frames} frames'
    # title += ', diff' if FRAME_DIFF else ', nodiff'
    # title += ', bkg rem' if REMOVE_BKG else ', no bkg rem'

    # a = np.abs(fourier)
    a = np.abs(scipy.fft.fftshift(fourier))**2
    # TODO: should we be doing fftshift in the DDM code?
    # print(a.mean(), a.std())
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    # ax.semilogv()
    log_a = np.log(a)
    # im = ax.imshow(log_a, vmin=log_a.mean()-2*log_a.std(), vmax=log_a.mean()+2*log_a.std(), extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
    if SAVE_TWO_D_PLOT:
        im = ax.imshow(log_a, extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
        fig.colorbar(im)
        ax.set_title(title)
        # ax.set_xlabel('$k_x$')
        # ax.set_ylabel('$k_y$')
        filename = f'static_fourier_{file}'
        common.save_fig(fig, f'DDM/figures_png/{filename}_old.png')

    ##################################
    # radial average
    ##################################
    title = common.name(file)
    title += f'\n{num_frames} frames'
    
    f = scipy.fft.fftshift( np.sqrt( fx**2 + fy**2 ) )
    print('f', f.min(), f.max())
    bins = np.linspace(0, f.max()/np.sqrt(2), 1000)[1:]
    print('bin_sep', 2*np.pi/(bins[1]-bins[0]))
    ff = np.abs(fourier)**2
    print('f, a', f.shape, a.shape)
    f_flat = f.flatten()
    a_flat = a.flatten()
    v, bins, _ = scipy.stats.binned_statistic(f_flat, a_flat, bins=bins)
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
    
    # ax.set_title(title + data.get('NAME', ''))
    # average_I_sq = image.mean()**2
    # print(average_I_sq, '<I^2>')

    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/static_fourier_av_{file}.pdf', hide_metadata=True)
    
    if particle_diameter:
        Pq = form_factor_sphere(k, particle_diameter/2)
        index = np.argmax(k > 1e0)
        Pq *= v[index] / Pq[index]
        ax.plot(k, Pq, label=f'$P(q)$ ($2R={particle_diameter:.0f}\mathrm{{\mu m}}$)')
        print(form_factor_sphere(k, 4))

    common.add_exponential_index_indicator(ax, exponent=-2, anchor=(2, 20), xlabel='k')
    
    plt.legend()

    filename = f'static_fourier_av_{file}'
    common.save_fig(fig, f'DDM/figures_png/{filename}_old.png', dpi=200)

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'static_fourier_'):

        show_static_fourier(file)