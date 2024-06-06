import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm

def speed_string(time_mult, every_nth_frame):
    if time_mult*every_nth_frame == 1:
        return 'realtime'
    return f'{time_mult*every_nth_frame}x speed'

def save_array_movie(stack, pixel_size, time_step, file, outputfilename, func=lambda timestep, ax : None):
    dpi = 200
    figsize = np.array(stack.shape)[[2, 1]] / dpi
    if figsize.mean() < 1.5:
        figsize *= 2
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    fig.tight_layout()
    fig.set_size_inches(*fig.get_size_inches()) # apparently this is needed to make subplots_adjust work
    ax.set_axis_off() # hide axes, ticks, ...
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

    time_mult = 4
    if file.startswith('marine'):
        time_mult = 0.25
    
    every_nth_frame = 1
    if file == 'pierre_exp':
        every_nth_frame = 10
        time_mult = 0.1

    fps = 1/time_step * time_mult
    frames = range(0, min(stack.shape[0], 50*every_nth_frame), every_nth_frame)


    def show(timestep):
        ax.clear()
        im = ax.imshow(stack[timestep, :, :], vmin=stack.min(), vmax=stack.max(), cmap=matplotlib.cm.Greys, interpolation='none')
        # if timestep == 0:
        #     fig.colorbar(im)
        # plt.imshow(stack.min(axis=0))
        ax.set_axis_off() # hide axes, ticks, ...
    
        common.add_scale_bar(ax, pixel_size)
        ax.text(0.95, 0.05, speed_string(time_mult, every_nth_frame), transform=ax.transAxes, ha='right')

        func(timestep, ax)
        # print(stack[:, :, timestep].mean())

    common.save_gif(show, frames, fig, outputfilename, fps=fps)

if __name__ == '__main__':
    for file in sys.argv[1:]:
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        stack = data['stack']
        pixel_size = data['pixel_size']
        time_step = data['time_step']

        # crop
        stack = stack[:, :500, :500]

        # stack = common.add_drift_intensity(stack, 1)

        print(stack.shape[1], 'x', stack.shape[2], 'px')

        for bkg in [True, False]:

            stack_copy = stack

            if bkg:
                stack_copy = stack_copy - stack_copy.mean(axis=0) # remove space background
            # print(stack.mean(axis=(1, 2)).std())
            # stack = stack - stack.mean(axis=(1, 2))[:, np.newaxis, np.newaxis] # remove total intensity fluctuations in time
                
            # print(stack.mean(axis=(1, 2)).std())

            filename = f'stack_movie_bkgrem_{file}' if bkg else f'stack_movie_{file}'
            save_array_movie(stack_copy, pixel_size, time_step, file, f"preprocessing/figures_png/{filename}.gif")
            # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4")