import sys
import matplotlib.pyplot as plt
import common
from scattering_functions.show_both import show_single_F_type

num_displayed_ks = 20
fig, axes = plt.subplots(4, num_displayed_ks, figsize=(num_displayed_ks*3, 4*2.8))

for i, file in enumerate(sys.argv[1:]):

    show_single_F_type(file, i, 'f', fig, axes, num_displayed_ks, mult=True)
            
    # plt.suptitle(fr'f(k, t), {file}')

filestring = '_'.join(sys.argv[1:])
common.save_fig(fig, f'scattering_functions/figures_png/f_mult_{filestring}.png', dpi=300)