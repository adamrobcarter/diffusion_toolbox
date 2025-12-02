import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft

FIRST_FRAME = False
FRAME_DIFF = False
REMOVE_BKG = False

if __name__ == '__main__':
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