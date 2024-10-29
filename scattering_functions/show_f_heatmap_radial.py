import common
import matplotlib.pyplot as plt
import numpy as np
import countoscope_theory.structure_factor
import scipy.optimize, scipy.stats
import functools
import matplotlib.cm

# SHOW_R_AXIS = False
# FORCE_LOG_X_AXIS = True
# RESCALE_X_AXIS_BY_DIAMETER = False
# SPLIT_AND_COLOR = False
# SHOW_FIT = False
# SHOW_THEORY = True
    
def show_at_t(time, ax, k_x, k_y, f, t, vmin, vmax, ylim, xlim):
    ax.clear()
    # ax.pcolormesh(k_x, k_y, np.log10(f[time, :, :]), shading='nearest', vmin=vmin, vmax=vmax)
    # im = ax.pcolormesh(k_x, k_y, f[time, :, :], shading='nearest', vmin=vmin, vmax=vmax)

    # ax.scatter(k_x.flatten(), k_y.flatten(), c=f[time, :, :].flatten())

    theta = np.arctan2(k_x, k_y) + np.pi/2
    r = np.sqrt(k_y**2 + k_x**2)
    # print(theta.max(), theta.min())
    assert np.all(r >= 0)

    im = ax.pcolormesh(theta, r, f[time, :, :], shading='nearest', vmin=vmin, vmax=vmax)
    
    # for x in range(f.shape[1]):
    #     for y in range(f.shape[2]):
    #         print(matplotlib.cm.viridis(f[time, x, y]))
    #         ax.add_patch(plt.Circle((k_x[x, y], k_y[x, y]), 0.1, color=matplotlib.cm.viridis(f[time, x, y])))

    ax.text(0.01, 0.01, f't={t[time]:.1f}', transform=ax.transAxes)
    ax.set_aspect('equal')
    ax.grid(False)
    # ax.set_xlim(*xlim)
    # ax.set_ylim(*ylim)
    # ax.set_xlabel('$L_x$')
    # ax.set_ylabel('$L_y$')
    # ax.set_xlabel('$k_x$')
    # ax.set_ylabel('$k_y$')

    # ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'${2*np.pi/x/particle_diameter:0.1f}\sigma$' if x != 0 else ''))
    # ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'${2*np.pi/x/particle_diameter:0.1f}\sigma$' if x != 0 else ''))
    # ax.set_xticks([-2, -1, -0.5, 0.5, 1, 2])

    return im# for colorbar etc

def go(file, export_destination=None):
    data = common.load(f"scattering_functions/data/F_{file}.npz")
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


    # # we specify the edges of the boxes
    # print(k_x.shape, k_y.shape, f.shape)
    # spacing_x = k_x[-1] - k_x[-2]
    # k_x = np.concatenate([k_x, [k_x[-1]+spacing_x]])
    # spacing_y = k_y[-1] - k_y[-2]
    # k_y = np.concatenate([k_y, [k_y[-1]+spacing_y]])
    print(k_x.shape, k_y.shape, f.shape)


    # need x and y to be a mesh
    # k_x_size = k_x.size
    # k_y_size = k_y.size
    # k_x = np.repeat(k_x[:, np.newaxis], k_y_size, axis=1)
    # k_y = np.repeat(k_y[np.newaxis, :], k_x_size, axis=0)

    # x = 2*np.pi/k_x / particle_diameter
    # y = 2*np.pi/k_y / particle_diameter
    # print(x)

    
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

    func = functools.partial(show_at_t, ax=ax, k_x=k_x, k_y=k_y, f=f, t=t, vmin=vmin, vmax=vmax, xlim=xlim, ylim=ylim)

    skip_frames = 1
    common.save_gif(func, range(0, t.size//skip_frames), fig, f'scattering_functions/figures_png/f_heatmap_{file}.gif', fps=4/skip_frames)

    # if export_destination:
    #     common.save_fig(fig, export_destination, hide_metadata=True)
    # common.save_fig(fig, f'scattering_functions/figures_png/S_of_k_all_{file}.png')



if __name__ == '__main__':
    for file in common.files_from_argv('scattering_functions/data', 'F_'):
        go(file)