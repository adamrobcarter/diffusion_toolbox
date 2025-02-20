import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

SAVE_TWO_D_PLOT = True
NUM_BINS = 50
NUM_KS = 10

def form_factor_sphere(q, R):
    Fs = 3 * (np.sin(q*R) - q*R*np.cos(q*R)) / (q*R)**3
    return Fs**2

def show_static_fourier(file):

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
    k_mag = np.sqrt(k_x**2 + k_y**2)
    k_ang = np.arctan2(k_x, k_y)

    title = common.name(file)
    title += f', {num_frames} frames'
    # title += ', diff' if FRAME_DIFF else ', nodiff'
    # title += ', bkg rem' if REMOVE_BKG else ', no bkg rem'

    # a = np.abs(fourier)
    fourier_abs = np.abs(fourier)**2
    # TODO: should we be doing fftshift in the DDM code?

    dk = 0.1

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    k_i = 0

    for k in np.logspace(np.log10(k_x[0, 1] - k_x[0, 0]), np.log10(np.nanmax(k_x)), num=NUM_KS):
        assert np.isfinite(k)
        mask = ( (k - dk/2) < k_mag ) & ( k_mag < (k + dk/2) )

        color = common.colormap(k_i, 0, NUM_KS)

        value, _, _ = scipy.stats.binned_statistic(k_ang[mask], fourier_abs[mask], statistic='mean', bins=NUM_BINS)
        angle, _, _ = scipy.stats.binned_statistic(k_ang[mask], k_ang[mask], statistic='mean', bins=NUM_BINS)

        um = r'\mathrm{\mu m^{-1}}'
        plt.plot(angle, value, label=f'$k={k:.3g}{um}$', color=color)

        k_i += 1

    ax.legend(fontsize=8)
    ax.semilogy()

    ax.set_xlabel('$\phi$')
    ax.set_ylabel(r'$\langle I(k, \phi)^* I(k, \phi) \rangle$')

    filename = f'static_fourier_phi_{file}'
    common.save_fig(fig, f'DDM/figures_png/{filename}.png', dpi=200)

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'static_fourier_'):

        show_static_fourier(file)