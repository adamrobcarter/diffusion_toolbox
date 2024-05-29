import common
import matplotlib.pyplot as plt
import numpy as np

for file in common.files_from_argv('scattering_functions/data', 'f_'):
    data = common.load(f"scattering_functions/data/F_{file}.npz")
    t                 = data["t"]
    F                 = data["F"] # (num timesteps) x (num k bins)
    F_unc             = data['F_unc']
    k                 = data["k"]
    particle_diameter = data['particle_diameter']

    S = F[0, :]

    fig, ax = plt.subplots(1, 1, figsize=(3.2, 3))
    
    start_index = 40
    ax.errorbar(k[0, start_index:], S[start_index:], yerr=F_unc[0, start_index:], linestyle='none', marker='o')
    # ax.loglog()
    # ax.semilogy()
    ax.set_ylim(0, 2)

    ax.set_ylabel('$S(k)$')
    ax.set_xlabel('$k$')

    
    if particle_diameter := data.get('particle_diameter'):
        # ax.vlines(2*np.pi/particle_diameter, S.min(), S.max())
        pass
    # ax.set_ylim(v.min()/1.1, v.max()*1.1)

    realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
    realspace_ax.set_xticks([1e1, 2e0, 1e0, 5e-1])
    realspace_ax.set_xlabel(r'$2\pi/k$ ($\mathrm{\mu m}$)')


    ax.text(0.7, 0.1, f'$\phi={file[-4:]}$', transform=ax.transAxes)

    common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/S_of_k_{file}.pdf', hide_metadata=True)
    common.save_fig(fig, f'scattering_functions/figures_png/S_of_k_{file}.png')


