import sys
import matplotlib.pyplot as plt
import common
from scattering_functions.show_both import show_single_F_type

for file in sys.argv[1:]:

    show_single_F_type(file, 0, 'f')
            
    plt.suptitle(fr'f(k, t), {file}')

    common.save_fig(plt.gcf(), f'scattering_functions/figures_png/Fs_decay_t_{file}.png', dpi=300)
        