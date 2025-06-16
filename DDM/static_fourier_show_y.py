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

def show_static_fourier(file, ax, color=None):

    data = common.load(f'DDM/data/static_fourier_{file}.npz')
    pixel_size        = data['pixel_size']
    num_frames        = data['num_frames_used']
    particle_diameter = data.get('particle_diameter')
    fourier_mod_sq    = data['fourier_mod_sq']
    k_x                = data['fx'] * 2 * np.pi
    k_y                = data['fy'] * 2 * np.pi

    fourier_mod_sq = np.fft.fftshift(fourier_mod_sq)
    k_x = np.fft.fftshift(k_x)
    k_y = np.fft.fftshift(k_y)

    title = common.name(file)
    title += f', {num_frames} frames'
    # title += ', diff' if FRAME_DIFF else ', nodiff'
    # title += ', bkg rem' if REMOVE_BKG else ', no bkg rem'

    # a = np.abs(fourier)
    fourier_abs = fourier_mod_sq
    # TODO: should we be doing fftshift in the DDM code?
    # print(a.mean(), a.std())
    # ax.semilogv()
    log_fourier_abs = np.log(fourier_abs)
    # log_fourier_abs = fourier_abs
    # im = ax.imshow(log_a, vmin=log_a.mean()-2*log_a.std(), vmax=log_a.mean()+2*log_a.std(), extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')



    k_x_is_zero = k_x == 0
    fourier_slice = fourier_abs[k_x_is_zero]
    k_y_slice = k_y[k_x_is_zero]
    print(fourier_slice.shape)


    ax.scatter(k_y_slice, fourier_slice, marker='.', color=color, s=5, label=file)
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel(r'$k$ ($\mathrm{\mu m}^{-1}$)')
    # ax.set_ylabel(r'$\langle I(k)^* I(k) \rangle$')

    # if particle_diameter:
    #     print(particle_diameter)
    #     print('onit', 2*np.pi/particle_diameter)
    #     ax.vlines(2*np.pi/particle_diameter, v.min(), v.max())
    # ax.vlines(2*np.pi/data['pixel_size'], v.min(), v.max())


    # ax.grid(alpha=0.3, which='both', axis='y')
    ax.grid(alpha=0.3)

    # end = 2e-1
    # end_index = np.argmax(k > end)
    # print(end_index)
    # func = lambda x, a, b, c: x**a
    # popt, pcov = scipy.optimize.curve_fit(func, k[:end_index], v[:end_index])
    # plt.plot(k[:end_index], func(k[:end_index], *popt), color='black', label=f'k**{popt[0]:.3f}')
    
    # func = lambda x, b, a: b * x + a
    # popt, pcov = scipy.optimize.curve_fit(func, np.log10(k[:end_index]), np.log10(v[:end_index]))
    # plt.plot(k[:end_index], 10**func(np.log10(k[:end_index]), *popt), color='black', label=f'k**{popt[0]:.3f}')
    

    realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
    # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
    realspace_ax.set_xlabel(r'$2\pi/k$ ($\mathrm{\mu m}$)')
    

def go(file):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    show_static_fourier(file, ax)

    filename = f'static_fourier_y_{file}'
    common.save_fig(fig, f'DDM/figures_png/{filename}.png', dpi=200)

def go_mult(files):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i, file in enumerate(files):
        show_static_fourier(file, ax, color=common.tab_color(i))

    filenames = '_'.join(files)
    filename = f'static_fourier_y_{filenames}'
    common.save_fig(fig, f'DDM/figures_png/{filename}.png', dpi=200)

if __name__ == '__main__':
    # for file in common.files_from_argv('DDM/data', 'static_fourier_'):
    #     go(file, SAVE_TWO_D_PLOT)

    go_mult(common.files_from_argv('DDM/data', 'static_fourier_'))