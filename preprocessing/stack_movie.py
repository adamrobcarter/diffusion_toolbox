import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm
import warnings

# HIGHLIGHTS = True # displays 50 frames evenly throughout the stack instead of the first 50
HIGHLIGHTS = False
# BACKWARDS = True
BACKWARDS = False # plays the stack backwards. DIFF_WITH_ZERO is now compared to the last frame

# possible display methods
REMOVE_BACKGROUND = 1
DIFF_WITH_ZERO = 2 # frame - frame1
DIFF_WITH_PREVIOUS = 3
NONE = 4

# select your method here
# METHOD = DIFF_WITH_PREVIOUS
# METHOD = DIFF_WITH_ZERO
METHOD = NONE
# METHOD = REMOVE_BACKGROUND

if METHOD == NONE:
    BACKWARDS = False

def speed_string(time_mult, every_nth_frame):
    if time_mult*every_nth_frame == 1:
        return 'realtime'
    return f'{time_mult*every_nth_frame}x speed'

def save_array_movie(stack, pixel_size, time_step, file, outputfilename, func=lambda timestep, ax : None, nth_frame=1):
    print('arrived in save_array_movie')
    dpi = 200
    figsize = np.array(stack.shape)[[2, 1]] / dpi
    while figsize.mean() < 3:
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
    print(fps, 'fps')

    # these are cause of the reduce operation
    if file.endswith('_25'):
        fps *= 25
    fps *= nth_frame#*every_nth_frame

    print(fps, 'fps', nth_frame)


    if fps < 0.2:
        warnings.warn(f'fps = {fps}')
    frames = range(0, min(stack.shape[0], 50*every_nth_frame), every_nth_frame)


    if BACKWARDS:
        stack = stack[::-1, :, :]

    if METHOD == REMOVE_BACKGROUND:
        print('subtracting mean')
        usedstack = stack[frames, :, :] - stack[:, :, :].mean(axis=0) # remove space background

    elif METHOD == DIFF_WITH_ZERO:
        print('subtracting frame 0')
        usedstack = stack[frames, :, :] - stack[0, :, :]

    elif METHOD == DIFF_WITH_PREVIOUS:
        print('subtracting previous frame')
        usedstack = stack[np.array(frames[:-1])+1, :, :] - stack[frames[:-1], :, :]

    if METHOD == NONE:
        usedstack = stack[frames, :, :]

    print('finding quantiles')
    vmin = np.quantile(usedstack, 0.01)
    vmax = np.quantile(usedstack, 0.99)
    # vmin = stack.min()
    # vmax = stack.max()
    # vmin = -5.5e-2
    # vmax = 6.2e-2
    # vmin = None
    # vmax = None
    print('ready to render')

    common.term_hist(usedstack)

    # print('min mean max', usedstack.min(), usedstack.mean(), usedstack.max())

    def show(index):
        if BACKWARDS:
            timestep = frames[-index]
        else:
            timestep = frames[index]
        ax.clear()
        im = ax.imshow(usedstack[index, :, :], vmin=vmin, vmax=vmax, cmap=matplotlib.cm.Greys_r, interpolation='none')
        # if timestep == 0:
        #     fig.colorbar(im)
        # plt.imshow(stack.min(axis=0))
        ax.set_axis_off() # hide axes, ticks, ...
    
    
        color = 'white' if usedstack.mean()/(usedstack.max()-usedstack.min()) < 0.2 else 'black'
        common.add_scale_bar(ax, pixel_size, color=color)
        ax.text(0.95, 0.1, speed_string(time_mult, every_nth_frame*nth_frame)+f'\ntime = {timestep*time_step*nth_frame:.1f}s', color=color, transform=ax.transAxes, ha='right', fontsize=15)
        
        func(timestep, ax)
        # print(stack[:, :, timestep].mean())

    common.save_gif(show, range(usedstack.shape[0]), fig, outputfilename, fps=fps)

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        stack = data['stack']
        pixel_size = data['pixel_size']
        time_step = data['time_step']

        # crop
        # stack = stack[:, :500, :500]

        # stack = common.add_drift_intensity(stack, 1)

        print(stack.shape[1], 'x', stack.shape[2], 'px')
        # print(stack.min(), stack.max())

        filename = f'stack_movie_{file}'
        if METHOD == REMOVE_BACKGROUND:
            filename += '_bkgrem'
        elif METHOD == DIFF_WITH_ZERO:
            filename += '_diffframe1'
        elif METHOD == DIFF_WITH_PREVIOUS:
            filename += '_diffprev'
        if HIGHLIGHTS:
            filename += '_highlights'
        if BACKWARDS:
            filename += '_backwards'
        save_array_movie(stack, pixel_size, time_step, file, f"preprocessing/figures_png/{filename}.gif",
                         nth_frame=data.get('nth_frame', 1))
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4")