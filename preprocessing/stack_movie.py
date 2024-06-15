import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm

HIGHLIGHTS = True

def speed_string(time_mult, every_nth_frame):
    if time_mult*every_nth_frame == 1:
        return 'realtime'
    return f'{time_mult*every_nth_frame}x speed'

def save_array_movie(stack, pixel_size, time_step, file, outputfilename, func=lambda timestep, ax : None, remove_background=True):
    print('arrived in save_array_movie')
    dpi = 200
    figsize = np.array(stack.shape)[[2, 1]] / dpi
    if figsize.mean() < 1.5:
        figsize *= 2
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    fig.tight_layout()
    fig.set_size_inches(*fig.get_size_inches()) # apparently this is needed to make subplots_adjust work
    ax.set_axis_off() # hide axes, ticks, ...
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

    print('set up figures')

    time_mult = 4
    if file.startswith('marine'):
        time_mult = 0.25
    
    if HIGHLIGHTS:
        every_nth_frame = stack.shape[0] // 50
    else:
        every_nth_frame = 1

    if file == 'pierre_exp':
        every_nth_frame = 10
        time_mult = 0.1

    fps = 1/time_step * time_mult
    frames = range(0, min(stack.shape[0], 50*every_nth_frame), every_nth_frame)

    print('finding quantiles')
    vmin = np.quantile(stack[frames, :, :], 0.05)
    vmax = np.quantile(stack[frames, :, :], 0.95)
    # vmin = stack.min()
    # vmax = stack.max()
    print('ready to render')

    

    if remove_background:
        print('subtracting mean')
        stack = stack - stack[:, :, :].mean(axis=0) # remove space background

    def show(timestep):
        ax.clear()
        im = ax.imshow(stack[timestep, :, :], vmin=vmin, vmax=vmax, cmap=matplotlib.cm.Greys, interpolation='none')
        # if timestep == 0:
        #     fig.colorbar(im)
        # plt.imshow(stack.min(axis=0))
        ax.set_axis_off() # hide axes, ticks, ...
    
        common.add_scale_bar(ax, pixel_size)
        ax.text(0.95, 0.05, speed_string(time_mult, every_nth_frame), transform=ax.transAxes, ha='right', fontsize=15)
        ax.text(0.95, 0.10, f'time = {int(timestep*time_step)} s',    transform=ax.transAxes, ha='right', fontsize=15)

        func(timestep, ax)
        # print(stack[:, :, timestep].mean())

    common.save_gif(show, frames, fig, outputfilename, fps=fps)

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        print(file)
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        stack = data['stack']
        pixel_size = data['pixel_size']
        time_step = data['time_step']

        # crop
        # stack = stack[:, :500, :500]

        # stack = common.add_drift_intensity(stack, 1)

        print(stack.shape[1], 'x', stack.shape[2], 'px')
        # print(stack.min(), stack.max())

        remove_background = True

        filename = f'stack_movie_{file}_bkgrem' if remove_background else f'stack_movie_{file}'
        if HIGHLIGHTS:
            filename += '_highlights'
        save_array_movie(stack, pixel_size, time_step, file, f"preprocessing/figures_png/{filename}.gif", remove_background=remove_background)
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4")