import common
import matplotlib.pyplot as plt
import numpy as np
import countoscope_theory.structure_factor
import scipy.optimize, scipy.stats
import isf.show_f_heatmap

TIME = 18
DO_F_NOT_f = False

# SMALL = False
# SHOW_R_AXIS = False
# FORCE_LOG_X_AXIS = True
# RESCALE_X_AXIS_BY_DIAMETER = False
# SPLIT_AND_COLOR = False
# SHOW_FIT = False
# SHOW_THEORY = True

def go(file, export_destination=None):
    data = common.load(f"isf/data/F_{file}.npz")
    t                 = data["t"]
    F                 = data["F_unbinned"] # (num timesteps) x (num kx) x (num ky)
    # F_unc             = data['F_unc']
    F_unc = np.zeros_like(F)
    k                 = data["k_unbinned"]
    k_x               = data["k_x"]
    k_y               = data["k_y"]
    particle_diameter = data.get('particle_diameter')

    S     = F    [0, :]
    k     = k    [0, :]
    S_unc = F_unc[0, :]

    if DO_F_NOT_f:
        f = F
        func_type = 'F'
    else:
        f = F / S
        func_type = 'f'
    # f = F


    print('asd', k_x[(k_x.size)//2], k_y[(k_y.size)//2], f[1, (k_x.size)//2, (k_y.size)//2])


    # need x and y to be a mesh
    k_x_size = k_x.size
    k_y_size = k_y.size
    k_x = np.repeat(k_x[:, np.newaxis], k_y_size, axis=1)
    k_y = np.repeat(k_y[np.newaxis, :], k_x_size, axis=0)

    # vmin = -5
    # vmax = 0
    # vmin = 0
    # vmax = 0.5
    vmin = 0
    vmax = 1

    axis_frac = 0.03
    xlim = (axis_frac*k_x.min(), axis_frac*k_x.max())
    ylim = (axis_frac*k_y.min(), axis_frac*k_y.max())
    print(xlim, ylim, 'xy')

    figsize = (0.1*(k_x.max()-k_x.min()), 0.1*(k_y.max()-k_y.min()))
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if t.size == 2:
        time_index = 1
    else:
        time_index = TIME
    
    title = f'${func_type}(k, {t[time_index]:.1f}s)$'

    if 'rot45' in file:
        title += ' rotated 45deg'
    if 'rot90' in file:
        title += ' rotated 90deg'

    # plt.suptitle(title) # idk why ax.set_title isn't working

    im = isf.show_f_heatmap.show_at_t(time_index, ax, k_x, k_y, f, t, vmin, vmax, xlim=xlim, ylim=ylim)

    # fig.colorbar(im)

    common.save_fig(fig, f'isf/figures_png/f_heatmap_{file}_singleframe.png')

    # if export_destination:
    #     common.save_fig(fig, export_destination, hide_metadata=True)
    # common.save_fig(fig, f'isf/figures_png/S_of_k_all_{file}.png')



if __name__ == '__main__':
    for file in common.files_from_argv('isf/data', 'F_'):
        go(file)