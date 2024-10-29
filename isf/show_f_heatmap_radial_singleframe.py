import common
import matplotlib.pyplot as plt
import numpy as np
from scattering_functions.show_f_heatmap_radial import show_at_t

def go(file, export_destination=None):
    data = common.load(f"isf/data/F_{file}.npz")
    t                 = data["t"]
    F                 = data["F_unbinned"] # (num timesteps) x (num kx) x (num ky)
    # F_unc             = data['F_unc']
    F_unc = np.zeros_like(F)
    k_x               = data["k_x"]
    k_y               = data["k_y"]
    particle_diameter = data.get('particle_diameter')

    S     = F    [0, :]
    # k     = k    [0, :]
    S_unc = F_unc[0, :]

    f = F / S
    # f = F

    # vmin = -5
    # vmax = 0
    vmin = 0
    vmax = 1

    
    axis_frac = 1
    xlim = (axis_frac*k_x.min(), axis_frac*k_x.max())
    ylim = (axis_frac*k_y.min(), axis_frac*k_y.max())

    figratio = (k_x.max()-k_x.min()) / (k_y.max()-k_y.min())
    figsize = 4 * np.array([np.sqrt(figratio), 1/np.sqrt(figratio)])
    print('figsize', figsize)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=figsize)
    ax.set_title('$F(k, t)$')

    show_at_t(1, ax, k_x, k_y, f, t, vmin=vmin, vmax=vmax, xlim=xlim, ylim=ylim)

    if export_destination:
        common.save_fig(fig, export_destination, hide_metadata=True)
    common.save_fig(fig, f'isf/figures_png/f_heatmap_{file}.png')



if __name__ == '__main__':
    for file in common.files_from_argv('isf/data', 'F_'):
        go(file)