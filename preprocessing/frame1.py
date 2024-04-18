import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm

SMALL = True

for file in common.files_from_argv('preprocessing/data/', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    # data = common.load('data/alice_stack_0.02.npz')
    # data = common.load('data/stack_alice_0.66.npz')
    stack = data['stack']

    print(stack.shape[1], 'x', stack.shape[2], 'px')


        
    # stack = stack - stack.mean(axis=0) # remove space background
    # stack = stack.mean(axis=0)
    plt.figure(figsize=(2.3, 2.3))
    
    frame1 = stack[0, :, :]

    cmap = matplotlib.cm.Greys
    if file.startswith('marine'):
        cmap = matplotlib.cm.Greys_r

    plt.imshow(frame1, cmap=cmap)
    # plt.imshow(stack.min(axis=0))
    
    # excess = stack - stack.min(axis=0)
    # print(stack[:, :, timestep].mean())
    print(file, frame1.mean()/(frame1.max()-frame1.min()))
    color = 'white' if frame1.mean()/(frame1.max()-frame1.min()) < 0.3 else 'black'

    if SMALL and file.startswith('eleanor0.'):
        plt.ylim(000, 600)
        plt.xlim(000, 600)
    if SMALL and file == 'pierre_exp':
        plt.ylim(000, 300)
        plt.xlim(000, 300)
    if SMALL and file.startswith('marine'):
        plt.ylim(000, 150)
        plt.xlim(000, 150)

    data2 = common.load(f'particle_detection/data/particles_{file}.npz')
    sigma = data2.get('particle_diameter')
    if not sigma or np.isnan(sigma):
        sigma = data2['particle_diameter_calced']
        assert not np.isnan(sigma)
    label = fr'$\sigma={sigma:.2f}\mathrm{{\mu m}}$'
    plt.gca().text(0.55, 0.05, label, color=color, transform=plt.gca().transAxes,
                   horizontalalignment='left', verticalalignment='bottom')
    
    if pack_frac := data.get('pack_frac_given'):
        plt.gca().text(0.55, 0.15, f'$\phi={pack_frac:.2f}$', color=color, transform=plt.gca().transAxes,
                    horizontalalignment='left', verticalalignment='bottom')

    common.add_scale_bar(plt.gca(), data['pixel_size'], color=color)

    common.save_fig(plt.gcf(), f'preprocessing/figures_png/frame1_{file}.png', dpi=600, only_plot=True)
    common.save_fig(plt.gcf(), f'/home/acarter/presentations/intcha24/figures/frame1_{file}.png', dpi=600, only_plot=True)