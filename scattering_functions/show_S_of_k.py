import common
import matplotlib.pyplot as plt
import numpy as np

SMALL = False

for file in common.files_from_argv('scattering_functions/data', 'F_'):
    data = common.load(f"scattering_functions/data/F_{file}.npz")
    t                 = data["t"]
    F                 = data["F"] # (num timesteps) x (num k bins)
    F_unc             = data['F_unc']
    k                 = data["k"]
    particle_diameter = data.get('particle_diameter', 1)

    S = F[0, :]

    fig, ax = plt.subplots(1, 1, figsize=(3.2, 3) if SMALL else (4, 4))
    
    start_index = 0 # crops off k=0 delta fn
    end_index = 1000
    ax.errorbar(particle_diameter*k[0, start_index:end_index], S[start_index:end_index], yerr=F_unc[0, start_index:end_index], linestyle='none', marker='o')
    # ax.loglog()
    # ax.semilogy()
    ax.set_ylim(0, 2)

    ax.set_xlim(0.01, min(15, particle_diameter*k.max()))

    ax.set_ylabel('$S(k)$')
    ax.set_xlabel('$k\sigma$')

    realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k/particle_diameter, lambda r: 2*np.pi*particle_diameter/r))
    # realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 1/k/particle_diameter, lambda r: 1*particle_diameter/r))
    realspace_ax.set_xticks([1e1, 1e0, 5e-1])
    realspace_ax.set_xlabel(r'$r/\sigma = 2\pi/k\sigma$')


    ax.text(0.7, 0.1, f'$\phi={file[-4:]}$', transform=ax.transAxes)

    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/S_of_k_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'scattering_functions/figures_png/S_of_k_{file}.png')


