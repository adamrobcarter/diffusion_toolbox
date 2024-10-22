import sys
import matplotlib.pyplot as plt
import common
from scattering_functions.show_both import show_single_F_type
import numpy as np

num_displayed_ks = 5

def go(files):
    fig, axes = plt.subplots(4, num_displayed_ks+10, figsize=((num_displayed_ks+10)*3, 4*2.8))

    for i, file in enumerate(files):

        if 'timeorigins' in file:
            num_timeorigins = int(file.split('timeorigins')[1])
            color = common.colormap(np.log2(num_timeorigins), 0, 5)
        else:
            color = None

        show_single_F_type(
            file, i, 'f', fig, axes, num_displayed_ks,
            mult=True, do_fits=False, markersize=1, errorbar_alpha=0, bad_alpha=0.4, color=color
        )
                
        # plt.suptitle(fr'f(k, t), {file}')

    filestring = '_'.join(sys.argv[1:])
    common.save_fig(fig, f'scattering_functions/figures_png/f_mult_{filestring}.png', dpi=300)

if __name__ == '__main__':
    go(common.files_from_argv('scattering_functions/data', 'F_'))