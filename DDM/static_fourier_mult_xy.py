import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

FIRST_FRAME = False
FRAME_DIFF = False
REMOVE_BKG = False

# fig, ax = plt.subplots(1, 1, figsize=(5, 5))

# for file in (files := common.files_from_argv('preprocessing/data', 'stack_')):
#     data = common.load(f'preprocessing/data/stack_{file}.npz')
#     stack      = data['stack']
#     pixel_size = data['pixel_size']

#     particle_diameter = data.get('particle_diameter')

#     if REMOVE_BKG:
#         print('subtracting mean')
#         stack = stack - stack.mean(axis=0)

#     if FIRST_FRAME:
#         if FRAME_DIFF:
#             images = stack[[1], :, :] - stack[[0], :, :]
#         else:
#             images = stack[[0], :, :]
#     else:
#         if FRAME_DIFF:
#             images = stack[1::5, :, :] - stack[:-1:5, :, :]
#         else:
#             images = stack[::5, :, :]
#     del stack
#     num_frames = images.shape[0]
    
#     # image = image[::4, ::4]
#     # pixel_size *= 4

#     # rng = np.random.default_rng()
#     # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
#     #                     2 * np.pi * np.arange(200) / 34)
#     # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

#     print(images.shape)
#     fx, fy, fouriers = common.fourier_2D(images, pixel_size, (1, 2))
#     del images
#     # print(fouriers.shape)
#     fourier = fouriers.mean(axis=0)
#     # fx = fx.mean(axis=0)
#     # fy = fy.mean(axis=0)
#     print('f', fourier.shape, fx.shape)



#     # a = np.abs(fourier)
#     a = np.abs(scipy.fft.fftshift(fourier))**2
#     # TODO: should we be doing fftshift in the DDM code?
#     # print(a.mean(), a.std())

    
#     # radial average
#     f = scipy.fft.fftshift( np.sqrt( fx**2 + fy**2 ) )
#     print('f', f.min(), f.max())
#     bins = np.linspace(0, f.max()/np.sqrt(2), 1000)[1:]
#     print('bin_sep', 2*np.pi/(bins[1]-bins[0]))
#     ff = np.abs(fourier)**2
#     print('f, a', f.shape, a.shape)
#     f_flat = f.flatten()
#     a_flat = a.flatten()
#     v, bins, _ = scipy.stats.binned_statistic(f_flat, a_flat, bins=bins)
#     err, _, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), statistic='std', bins=bins)
#     n, _, _ = scipy.stats.binned_statistic(f.flatten(), a.flatten(), statistic='count', bins=bins)
#     vmean = np.mean(v[-20:-1])
#     x = (bins[1:] + bins[:-1]) / 2
#     k = 2 * np.pi * x
#     ax.scatter(k, v/vmean, s=2)
#     ax.errorbar(k, v/vmean, yerr=err/np.sqrt(n)/vmean, linestyle='none', label=common.name(file))

#     # ax.vlines(2*np.pi/1.5, *ax.get_ylim())

#     if particle_diameter:
#         print(particle_diameter)
#         print('onit', 2*np.pi/particle_diameter)
#         ax.vlines(2*np.pi/particle_diameter, v.min(), v.max())
#     # ax.vlines(2*np.pi/data['pixel_size'], v.min(), v.max())
#     # ax.set_ylim(v.min()/1.1, v.max()*1.1)

# ax.grid()
# ax.semilogy()
# ax.semilogx()
# ax.set_xlabel(r'$k$ ($\mathrm{\mu m}^{-1}$)')
# ax.set_ylabel(r'$\langle I(k)^* I(k) \rangle$')

# end = 2e-1
# end_index = np.argmax(k > end)
# print(end_index)
#     # func = lambda x, a, b, c: x**a
#     # popt, pcov = scipy.optimize.curve_fit(func, k[:end_index], v[:end_index])
#     # plt.plot(k[:end_index], func(k[:end_index], *popt), color='black', label=f'k**{popt[0]:.3f}')
    
#     # func = lambda x, b, a: b * x + a
#     # popt, pcov = scipy.optimize.curve_fit(func, np.log10(k[:end_index]), np.log10(v[:end_index]))
#     # plt.plot(k[:end_index], 10**func(np.log10(k[:end_index]), *popt), color='black', label=f'k**{popt[0]:.3f}')
    


# realspace_ax = ax.secondary_xaxis('top', functions=(lambda k: 2*np.pi/k, lambda r: 2*np.pi/r))
# # realspace_ax.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
# realspace_ax.set_xlabel(r'$2\pi/k$ ($\mathrm{\mu m}$)')

# ax.set_xlim(k.min()/1.1, max(k.max(), 2*np.pi/pixel_size)*1.1)

# ax.set_title(title)

# ax.legend(fontsize=7)
# # average_I_sq = image.mean()**2
# # print(average_I_sq, '<I^2>')

# # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/static_fourier_av_{file}.pdf', hide_metadata=True)

# filestr = '_'.join(files)
# filename = f'static_fourier_mult_{filestr}'
# filename += '_diff' if FRAME_DIFF else '_nodiff'
# filename += '_bkgrem' if REMOVE_BKG else '_nobkgrem'
# common.save_fig(fig, f'DDM/figures_png/{filename}.png', dpi=200)

##########################################
# now plot the x and y components respectively 
##########################################

cmap = plt.get_cmap('viridis')
iplot = 1
fig, ax = plt.subplots(1, 1, figsize=(5, 5))

for file in (files := common.files_from_argv('preprocessing/data', 'stack_')):
    
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']

    particle_diameter = data.get('particle_diameter')

    if REMOVE_BKG:
        print('subtracting mean')
        stack = stack - stack.mean(axis=0)

    if FIRST_FRAME:
        if FRAME_DIFF:
            images = stack[[1], :, :] - stack[[0], :, :]
        else:
            images = stack[[0], :, :]
    else:
        if FRAME_DIFF:
            images = stack[1::5, :, :] - stack[:-1:5, :, :]
        else:
            images = stack[::5, :, :]
    del stack
    num_frames = images.shape[0]
    
    # image = image[::4, ::4]
    # pixel_size *= 4

    # rng = np.random.default_rng()
    # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
    #                     2 * np.pi * np.arange(200) / 34)
    # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

    print(images.shape)
    fx, fy, fouriers = common.fourier_2D(images, pixel_size, (1, 2))
    del images
    # print(fouriers.shape)
    fourier = fouriers.mean(axis=0)
    # fx = fx.mean(axis=0)
    # fy = fy.mean(axis=0)
    print('f', fourier.shape, fx.shape)



    # a = np.abs(fourier)
    a = np.abs(scipy.fft.fftshift(fourier))**2
    log_a = np.log(a)

    lf = log_a.shape[1]
    aver = 5
    axt = np.zeros((aver,lf))
    ayt = np.zeros((aver,lf))
    for ix in range(aver): #this is stupid but I couldn't find a better way of doing this
         for iy in range(lf):
            axt[ix][iy] = log_a[lf//2 - aver//2 + ix][iy]
            ayt[ix][iy] = log_a[iy][lf//2 - aver//2 + ix]

    axx = axt.mean(axis = 0)
    ayy = ayt.mean(axis = 0)
    axx = axx[0:lf//2-1]
    ayy = ayy[0:lf//2-1]

    
    bins = np.linspace(fx.max(), 0, lf//2+1)[1:]
    x = (bins[1:] + bins[:-1]) / 2
    k = 2 * np.pi * x
    k.transpose()

    meanampx = np.mean(axx[-20:-1])
    meanampy = np.mean(ayy[-20:-1])
    ax.plot(k, axx/meanampx,'.', ms = 2, label = common.name(file)+'_x',color = cmap(iplot/len(files)))
    ax.plot(k, ayy/meanampy,'d', ms = 2,color = [c*0.9 for c in cmap(iplot/len(files))])
    iplot += 1

ax.semilogx()
ax.set_xlabel(r'$k$ ($\mathrm{\mu m}^{-1}$)')
ax.set_ylabel(r'$\langle I(kx,ky) \rangle$')

ax.grid()
ax.legend()


title = ''
title += f'\n{num_frames} frames'
title += ', diff' if FRAME_DIFF else ', nodiff'
title += ', bkg rem' if REMOVE_BKG else ', no bkg rem'
ax.set_title(title)
    
filestr = '_'.join(files)
filename = f'static_fourier_mult_xy_{filestr}'
filename += '_diff' if FRAME_DIFF else '_nodiff'
filename += '_bkgrem' if REMOVE_BKG else '_nobkgrem'
common.save_fig(fig, f'DDM/figures_png/{filename}_xy.png', dpi=200)