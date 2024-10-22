import sys
import matplotlib.pyplot as plt
import common

for file in common.files_from_argv('scattering_functions/data/', 'F_'):
    from scattering_functions.show_both import show_single_F_type
    # ^^^ rn we do this hack cause the figure is created at import inside show_both

    num_displayed_ks = 10
    fig, axes = plt.subplots(4, num_displayed_ks+5, figsize=((num_displayed_ks+5)*3, 4*2.8))
                                    # hack      ^^
    show_single_F_type(file, 0, 'f', fig, axes, num_displayed_ks)
            
    fig.suptitle(fr'f(k, t), {file}')

    common.save_fig(fig, f'scattering_functions/figures_png/f_decay_t_{file}.png', dpi=300)