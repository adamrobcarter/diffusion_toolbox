import common
import matplotlib.pyplot as plt
import numpy as np
import countoscope_theory.structure_factor
import scipy.optimize, scipy.stats
import functools

# SHOW_R_AXIS = False
# FORCE_LOG_X_AXIS = True
# RESCALE_X_AXIS_BY_DIAMETER = False
# SPLIT_AND_COLOR = False
# SHOW_FIT = False
# SHOW_THEORY = True
    
def show_at_t(time, fig, k_x, k_y, f, t, vmin, vmax, ylim, xlim):
    fig.clear()
    ax = fig.add_subplot() # the old axes gets deleted by fig.clear()
    # im = ax.pcolormesh(k_x, k_y, f[time, :, :], shading='nearest', vmin=vmin, vmax=vmax)

    data = np.fft.fftshift(f[time, :, :])
    k_x = np.fft.fftshift(k_x)
    k_y = np.fft.fftshift(k_y)

    im = ax.pcolormesh(k_x, k_y, data, shading='nearest', vmin=vmin, vmax=vmax)
    # im = ax.imshow(data, interpolation='none', vmin=vmin, vmax=vmax)
    fig.colorbar(im)

    ax.text(0.05, 0.05, f't={t[time]:.1f}', transform=ax.transAxes)
    ax.set_aspect('equal')
    # ax.set_xlim(*xlim)
    # ax.set_ylim(*ylim)
    # ax.set_xlabel('$L_x$')
    # ax.set_ylabel('$L_y$')
    ax.set_xlabel('$k_x$')
    ax.set_ylabel('$k_y$')

    # one_px_x = np.abs(k_x[0, 1] - k_x[0, 0])
    # one_px_y = np.abs(k_y[1, 0] - k_y[0, 0])
    # assert one_px_x > 0
    # assert one_px_y > 0
    # ax.set_xscale('symlog', linthresh=10*one_px_x)
    # ax.set_yscale('symlog', linthresh=10*one_px_y)

    # ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'${2*np.pi/x/particle_diameter:0.1f}\sigma$' if x != 0 else ''))
    # ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'${2*np.pi/x/particle_diameter:0.1f}\sigma$' if x != 0 else ''))
    # ax.set_xticks([-2, -1, -0.5, 0.5, 1, 2])

def go(file, export_destination=None):
    data = common.load(f"DDM/data/ddm_{file}.npz")
    t                 = data["t"]
    F                 = data["F_D_sq_noradial"] # (num timesteps) x (num kx) x (num ky)
    # F_unc             = data['F_unc']
    # F_unc = np.zeros_like(F)
    # k                 = data["k_unbinned"]
    k_x               = data["k_x"]
    k_y               = data["k_y"]

    # S     = F    [0, :]
    # k     = k    [0, :]
    # S_unc = F_unc[0, :]

    # f = F / S
    # f = F


    # # we specify the edges of the boxes
    # print(k_x.shape, k_y.shape, f.shape)
    # spacing_x = k_x[-1] - k_x[-2]
    # k_x = np.concatenate([k_x, [k_x[-1]+spacing_x]])
    # spacing_y = k_y[-1] - k_y[-2]
    # k_y = np.concatenate([k_y, [k_y[-1]+spacing_y]])


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
    # vmin = np.nanmin(F)
    # vmax = np.nanmax(F)

    F = F[1:, :, :] # t=0 doesn't make sense
    t = t[1:]

    F = np.log10(F)
    # F = F / F[0, :, :]

    vmin = np.nanmin(F)
    vmax = np.nanmax(F)
    # print(F)
    vmin = np.nanquantile(F, 0.002)
    vmax = np.nanquantile(F, 0.998)
    print(vmin, vmax)
    if 'psiche' in file:
        vmin = np.log10(1000)
        vmax = np.log10(5000)
    # vmin = None
    # vmax = None
    # vmin = 0
    assert np.isfinite(F).all()
    assert np.isfinite(vmin)
    assert np.isfinite(vmax)

    
    axis_frac = 0.03
    # xlim = (axis_frac*k_x.min(), axis_frac*k_x.max())
    # ylim = (axis_frac*k_y.min(), axis_frac*k_y.max())

    # figsize = (0.1*(k_x.max()-k_x.min()), 0.1*(k_y.max()-k_y.min()))
    # print('figsize', figsize)
    fig = plt.figure()#, figsize=figsize)
    # ax.set_title('$F(k, t)$')

    func = functools.partial(show_at_t, fig=fig, k_x=k_x, k_y=k_y, f=F, t=t, vmin=vmin, vmax=vmax, xlim=None, ylim=None)

    skip_frames = 1
    common.save_gif(func, range(0, t.size//skip_frames), fig, f'DDM/figures_png/DDM_heatmap_{file}.gif', fps=4/skip_frames)

    # if export_destination:
    #     common.save_fig(fig, export_destination, hide_metadata=True)
    # common.save_fig(fig, f'isf/figures_png/S_of_k_all_{file}.png')



if __name__ == '__main__':
    for file in common.files_from_argv('DDM/data', 'ddm_'):
        go(file)