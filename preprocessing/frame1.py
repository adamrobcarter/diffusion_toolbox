import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm

SMALL = False
REMOVE_BACKGROUND = False

for file in common.files_from_argv('preprocessing/data/', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    # data = common.load('data/alice_stack_0.02.npz')
    # data = common.load('data/stack_alice_0.66.npz')
    stack = data['stack']

    print(stack.shape[1], 'x', stack.shape[2], 'px')
    print(stack.shape[1]*data['pixel_size'], 'x', stack.shape[2]*data['pixel_size'], 'um')

        
    # stack = stack - stack.mean(axis=0) # remove space background
    # stack = stack.mean(axis=0)
    fig, ax = plt.subplots(1, 1, figsize=(2.3, 2.3) if SMALL else None)
    
    frame1 = stack[0, :, :]

    if REMOVE_BACKGROUND:
        print('subtracting background!!!')
        frame1 -= stack.mean(axis=0)

    cmap = matplotlib.cm.Greys
    if file.startswith('marine'):
        cmap = matplotlib.cm.Greys_r

    print(frame1.min(), frame1.max(), frame1.mean())
    sorted = np.sort(frame1.flatten())
    ax.imshow(frame1, cmap=cmap, interpolation='none')
    # plt.imshow(stack.min(axis=0))
    
    # excess = stack - stack.min(axis=0)
    # print(stack[:, :, timestep].mean())
    print(file, frame1.mean()/(frame1.max()-frame1.min()))

    color = 'white' if frame1.mean()/(frame1.max()-frame1.min()) < 0.2 else 'black'

    if SMALL and file.startswith('eleanor0.'):
        ax.ylim(000, min(400, stack.shape[2]))
        ax.xlim(000, min(400, stack.shape[1]))
    if SMALL and file == 'pierre_exp':
        ax.ylim(000, 300)
        ax.xlim(000, 300)
    if SMALL and file.startswith('marine'):
        ax.ylim(000, 150)
        ax.xlim(000, 150)

    try:
        data2 = common.load(f'particle_detection/data/particles_{file}.npz')
        sigma = data2.get('particle_diameter')
        if not sigma or np.isnan(sigma):
            sigma = data2['particle_diameter_calced']
            assert not np.isnan(sigma)
        label = fr'$\sigma={sigma:.2f}\mathrm{{\mu m}}$'
        ax.text(0.45, 0.05, label, color=color, transform=ax.transAxes,
                    horizontalalignment='left', verticalalignment='bottom')
        
        if pack_frac := data.get('pack_frac_given'):
            ax.text(0.45, 0.15, f'$\phi={pack_frac:.2f}$', color=color, transform=ax.transAxes,
                        horizontalalignment='left', verticalalignment='bottom')
    except:
        print('failed to find particle diamter / pack frac for annotations')

    common.add_scale_bar(ax, data['pixel_size'], color=color)

    filename = f'frame1_{file}'
    if REMOVE_BACKGROUND:
        filename += '_bkgrem'
    common.save_fig(fig, f'preprocessing/figures_png/{filename}.png', dpi=600, only_plot=True)
    # common.save_fig(fig, f'/home/acarter/presentations/cin_first/figures/frame1_{file}.png', dpi=300, only_plot=True)