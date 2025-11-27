import common
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data/', 'A_B_of_q_'):
        data = common.load(f'DDM/data/A_B_of_q_{file}.npz')
        A = data['A']
        B = data['B']
        q = data['q']
        pack_frac = data['pack_frac_given']
        sigma = data['particle_diameter']

        data2 = common.load(f'particle_detection/data/particle_intensity_profile_{file}.npz')
        profile = data2['profile']
        profile_k = data2['k']
        i = np.interp(q, profile_k, profile)

        F_k_0 = common.structure_factor_2d_hard_spheres(q, pack_frac, sigma)

        fig, (ax_A, ax_B, ax_F, ax_i) = plt.subplots(4, 1, figsize=(3, 10))

        # ax_A.scatter(q, A)
        ax_A.scatter(q, A/F_k_0)
        # ax_A.scatter(q, A/i**2/F_k_0)
        ax_B.scatter(q, B)
        ax_F.scatter(q, F_k_0)
        ax_i.scatter(profile_k, profile)
        ax_i.scatter(q, i)
        ax_i.semilogy()
        common.save_fig(plt.gcf(), f'DDM/figures_png/A_B_of_q_{file}.png')