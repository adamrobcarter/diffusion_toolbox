import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize


def do_intensity_means(file, stack, pixel_size, time_step):
    
    imeanint = 20
    shapes = stack.shape
    imean = np.zeros(shapes[0]-imeanint)
    for it in range(len(imean)):
        sslab = stack[it:it+imeanint, :, :]
        imean[it] = sslab.mean()
        

    title = common.name(file)
    title += f' mean intensity'

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.plot( [it*time_step for it in range(len(imean))],imean, '.')

    ax.set_xlabel(r'$t$ (s)')
    ax.set_ylabel(r'Rolling mean image intensity $\langle I(x,y,t) \rangle_{[t,t+20],x,y}$')

    # plt.imshow(image)
    ax.set_title(title)
    

    # ax.set_xlabel('$k_x$')
    # ax.set_ylabel('$k_y$')
    common.save_fig(fig, f'DDM/figures_png/intensity_means_{file}.png')

    imeanint = 20
    shapes = stack.shape
    imeandiff = np.zeros(shapes[0]-imeanint)
    sdiffs = stack[1:,:,:] - stack[:-1,:,:]
    for it in range(len(imean)):
        imeandiff[it] = sdiffs[it:it+imeanint-1, :, :].mean()
        

    title = common.name(file)
    title += f' mean intensity diff'

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.plot( [it*time_step for it in range(len(imean))],imeandiff, '.')

    ax.set_xlabel(r'$t$ (s)')
    ax.set_ylabel(r'Rolling mean diff intensity $\langle I(x,y,t+1) - I(x,y,t) \rangle_{[t,t+20],x,y}$')

    # plt.imshow(image)
    ax.set_title(title)
    

    # ax.set_xlabel('$k_x$')
    # ax.set_ylabel('$k_y$')
    common.save_fig(fig, f'DDM/figures_png/intensity_means_diff_{file}.png')

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data['stack']
        pixel_size = data['pixel_size']
        time_step  = data['time_step']

    do_intensity_means(file, stack, pixel_size, time_step)