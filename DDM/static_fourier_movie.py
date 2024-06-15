import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
import scipy.fft
import names

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        

    # images = stack[0, :, :]
    def show_n(n):
        ax.clear()
        images = stack[n, :, :]
        # if stack.shape[0] > 1:
        #     print('subbing back')
        #     image = image - stack.mean(axis=0)
        
        # image = image[::4, ::4]
        # pixel_size *= 4

        # rng = np.random.default_rng()
        # [X, Y] = np.meshgrid(2 * np.pi * np.arange(200) / 12,
        #                     2 * np.pi * np.arange(200) / 34)
        # image = np.sin(X) + np.cos(Y) + rng.uniform(0, 1, X.shape)*0.5

        print(images.shape)
        fx, fy, fourier = common.fourier_2D(images, pixel_size, (1, 2))
        # print(fouriers.shape)
        # fourier = fouriers.mean(axis=0)
        # fx = fx.mean(axis=0)
        # fy = fy.mean(axis=0)
        print('f', fourier.shape, fx.shape)

        # a = np.abs(fourier)
        a = np.abs(scipy.fft.fftshift(fourier))**2
        # TODO: should we be doing fftshift in the DDM code?
        # print(a.mean(), a.std())
        # ax.semilogv()
        log_a = np.log(a)
        im = ax.imshow(log_a, vmin=log_a.mean()-2*log_a.std(), vmax=log_a.mean()+2*log_a.std(), extent=(fx.min(), fx.max(), fy.min(), fy.max()), interpolation='none')
        plt.colorbar(im)
        # plt.imshow(image)
        ax.set_title(names.name(file))
    
    
    common.save_gif(show_n, range(50), fig, f'DDM/figures_png/static_fourier_{file}.gif', fps=4)
    