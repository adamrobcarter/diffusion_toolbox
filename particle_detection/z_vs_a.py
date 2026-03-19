import common
import matplotlib.pyplot as plt
import numpy as np

Z_OVER_A = False

if __name__ == '__main__':
    fig, ax = plt.subplots()

    z_mean = []
    z_std = []
    a = []

    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        particles = data['particles']
        a.append(data['particle_diameter'] / 2)
        z = particles[:, 2]

        z_mean.append(z.mean())
        z_std.append(z.std())

    z_mean = np.array(z_mean)
    z_std = np.array(z_std)
    a = np.array(a)

    for i in range(len(a)):
        print(f'a={a[i]}: <z> = {common.format_val_and_unc(z_mean[i]/a[i], z_std[i]/a[i], latex=False, sigfigs=3)}a')
    
    x_th = np.linspace(0, max(a))
    func = lambda x, m, b: m*x + b
    popt, pcov = common.curve_fit(func, a, z_mean)

    if Z_OVER_A:
        ax.errorbar(a, z_mean/a, yerr=z_std/a, marker='o', linestyle='none', label='simulation')

        ax.plot(x_th, func(x_th, *popt)/x_th, label=fr'fit: $\langle z \rangle = {popt[0]:.3f} a + {popt[1]:.3f}$')

    else:
        ax.errorbar(a, z_mean, yerr=z_std, marker='o', linestyle='none', label='simulation')


        ax.plot(x_th, func(x_th, *popt), label=fr'fit: $\langle z \rangle = {popt[0]:.3f} a + {popt[1]:.3f}$')

    ax.legend()
    if not Z_OVER_A:
        ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    ax.set_xlabel('$a$')
    if Z_OVER_A:
        ax.set_ylabel('$z/a$')
    else:
        ax.set_ylabel('$z$')
    phi = data['phi']
    ax.set_title(fr'Mean particle height ($\phi={phi}$)')
    filename = f'particle_detection/figures_png/z_vs_a_{file}'
    filename += '_zoa' if Z_OVER_A else ''
    filename += '.png'
    common.save_fig(fig, filename)