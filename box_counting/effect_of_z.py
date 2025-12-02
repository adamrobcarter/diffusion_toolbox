import matplotlib.pyplot as plt
import common
import numpy as np

# for L in [0.1, 1, 10]:
#     for z_mult in np.logspace(-0.5, 0.5, 5):
#         # L = 1
#         Lz = L * z_mult
#         D = 1
#         t = np.logspace(-4, 2)
#         N = 1
#         N2_theory = common.N2_nointer_3D(t, D, N, L, L, Lz)
#         plt.plot(t, N2_theory, label=f'Lz/L = {z_mult:.2f}')

if __name__ == '__main__':
    for D, z_mult in ((1, 1), (4, 1), (1, 1/4)):
        L = 1
        Lz = L * z_mult
        t = np.logspace(-4, 2)
        N = 1
        N2_theory = common.N2_nointer_3D(t, D, N, L, L, Lz)
        plt.plot(t, N2_theory, label=f'Lz/L = {z_mult:.2f}, D={D}')

    plt.loglog()
    plt.legend()
    common.save_fig(plt.gcf(), 'box_counting/figures_png/effect_of_z.png')