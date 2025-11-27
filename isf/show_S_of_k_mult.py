import common
import matplotlib.pyplot as plt
import numpy as np
import countoscope_theory.structure_factor

SMALL = False
SHOW_R_AXIS = False
SHOW_THEORY = True

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize=(3.2, 3) if SMALL else (4, 4))

    pack_frac = None
    particle_diameter_global = None

    def go(file, color):
        global pack_frac, particle_diameter_global, k

        data = common.load(f"isf/data/F_first_{file}.npz")
        t                 = data["t"]
        F                 = data["F"] # (num timesteps) x (num k bins)
        F_unc             = data['F_unc']
        k                 = data["k"]
        particle_diameter = data.get('particle_diameter', 1)

        S = F[0, :]

        start_index = 0 # crops off k=0 delta fn
        start_index = 0
        end_index = 1000
        ax.errorbar(particle_diameter*k[start_index:end_index], S[start_index:end_index], yerr=F_unc[0, start_index:end_index], linestyle=':', marker='o', label=file, markersize=4, color=color)

        if SHOW_R_AXIS:
            realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k/particle_diameter, lambda r: 2*np.pi*particle_diameter/r))
            realspace_ax.set_xticks([1e1, 1e0, 5e-1])
            realspace_ax.set_xlabel(r'$r/\sigma = 2\pi/k\sigma$')

        if SHOW_THEORY:
            if pack_frac:
                assert data['pack_frac'] == pack_frac
                assert particle_diameter == particle_diameter_global
            else:
                pack_frac = data['pack_frac']
                particle_diameter_global = particle_diameter


    for file_i, file in enumerate(filenames := common.files_from_argv('isf/data', 'F_first_')):
        color = common.colormap(file_i, 0, len(filenames))
        go(file, color)


    if SHOW_THEORY:
        ax.plot(particle_diameter_global*k, countoscope_theory.structure_factor.hard_spheres_2d(k, pack_frac, particle_diameter_global), color='black')

    ax.set_ylim(0.3, 1.2)
    ax.legend(fontsize=8)
    ax.semilogx()
    ax.set_xlim(0.1, 40)
    ax.set_ylabel('$S(k)$')
    ax.set_xlabel('$k\sigma$')

    # ax.text(0.7, 0.1, f'$\phi={file[-4:]}$', transform=ax.transAxes)

    filenames_str = '_'.join(filenames)

    # if export_destination:
    #     common.save_fig(fig, export_destination, hide_metadata=True)
    common.save_fig(fig, f'isf/figures_png/S_of_k_{filenames_str}.png', dpi=300)

    # if __name__ == '__main__':