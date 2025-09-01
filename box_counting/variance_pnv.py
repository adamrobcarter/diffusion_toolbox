import common
import matplotlib.pyplot as plt
import countoscope_theory.nmsd
import countoscope_theory.structure_factor
import numpy as np

for file in common.files_from_argv('box_counting/data', 'pnv_'):
    data = common.load(f'box_counting/data/pnv_{file}.npz')

    N_mean = data['N_mean']
    N_var  = data['N_var']
    N1N2mN1N2  = data['N1N2mN1N2']
    L      = data['box_sizes_x']
    phi    = data['pack_frac']
    sigma  = data['particle_diameter']

    fig, ax = plt.subplots(1, 1)
    ax.scatter(L, N_mean, label=r'$\langle N \rangle$')
    ax.scatter(L, N_var, label=r'$\mathrm{Var}(N)$')
    ax.scatter(L, N_mean * common.S_k_zero(phi), label=r'$\langle N \rangle S(k=0)$')
    ax.scatter(L, -N1N2mN1N2)#, label=r'$\langle N_1N_2 \rangle - \langle N_1 \rangle \langle N_2 \rangle$')

    theory_Ls = np.logspace(np.log10(L.min()), np.log10(L.max()))
    theory_plats = np.full_like(theory_Ls, np.nan)
    for i, theory_L in enumerate(theory_Ls):
        theory_plats[i] = 0.5 * countoscope_theory.nmsd.plateau_inter_2d(theory_L**2*data['density'], theory_L, lambda k: countoscope_theory.structure_factor.hard_spheres_2d(k, phi, sigma))
    ax.plot(theory_Ls, theory_plats, label='theory (PRX eq 6)', color='black')


    ax.legend()
    ax.semilogy()

    ax.set_xlabel(r'$L$')
    ax.semilogx()

    ax.set_title(f'$\phi$ = {phi:.2f}')

    common.save_fig(fig, f'box_counting/figures_png/variance_pnv_{file}.png')
