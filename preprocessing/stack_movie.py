import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm
import warnings

DISPLAY_SMALL = False
INVERSE_COLORS = False

# HIGHLIGHTS = True # displays 50 frames evenly throughout the stack instead of the first 50
HIGHLIGHTS = False
# BACKWARDS = True
BACKWARDS = False # plays the stack backwards. DIFF_WITH_ZERO is now compared to the last frame

# possible display methods
REMOVE_BACKGROUND = 1
DIFF_WITH_ZERO = 2 # frame - frame1
DIFF_WITH_PREVIOUS = 3
NONE = 4
TWOCHANNEL = 5

# select your method here
# METHOD = DIFF_WITH_PREVIOUS
# METHOD = DIFF_WITH_ZERO
METHOD = NONE
# METHOD = REMOVE_BACKGROUND

ADD_DRIFT = False

if METHOD == NONE:
    BACKWARDS = False

def normalise_stack(stack, file):
    max = np.quantile(stack, 0.99)
    stack = stack / max
    stack[stack > 1] = 1
    print(file, max)
    return stack

def speed_string(time_mult, every_nth_frame):
    if time_mult*every_nth_frame == 1:
        return 'realtime'
    return f'{time_mult*every_nth_frame}x speed'

def save_array_movie(stack, pixel_size, time_step, file, outputfilename,
                     window_size_x, window_size_y,
                     func=lambda timestep, ax : None, nth_frame=1, max_num_frames=50,
                     dpi=300,
                     display_small=True, inverse_colors=False, highlights=False,
                     backwards=False, method=NONE, stacks=None, stackcolors=None, channel=None,
                     dataset_name=None,
                     num_timesteps_in_data=None, # how many timesteps are in the original data?. only needed if stack is None
                     ):

    figsize = (3, 3)

    no_stack = stack == None

    if not no_stack:
        if method == TWOCHANNEL:
            assert stack is None
            assert stacks is not None
            stack = stacks[0] # this copy is so that we can use the same checks in twochannel mode as we do in one channnel mode
        else:
            assert stack is not None
            assert stacks is None

    
        figsize = np.array(stack.shape)[[2, 1]] / dpi
        while figsize.mean() < 3:
            figsize *= 2
        # if display_small and not : # potentially remove this for marine
        #     figsize = (2.3, 2.3)
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    time_mult = 4
    if file.startswith('marine'):
        time_mult = 0.25
    if file.startswith('marine'):
        time_mult = 0.25
    
    if highlights:
        every_nth_frame = stack.shape[0] // 50
    else:
        every_nth_frame = 1

    if file == 'pierre_exp':
        every_nth_frame = 10
        time_mult = 0.1

    fps = 1/time_step * time_mult

    # these are cause of the reduce operation
    if file.endswith('_25'):
        fps *= 25
    fps *= nth_frame#*every_nth_frame


    if fps < 0.2:
        warnings.warn(f'fps = {fps}')

    if not no_stack:
        num_timesteps_in_data = stack.shape[0]
    frames = range(0, min(num_timesteps_in_data, max_num_frames*every_nth_frame), every_nth_frame)

    if not no_stack:
        if backwards:
            raise # remove me whenever you like
            stack = stack[::-1, :, :]

        if method == REMOVE_BACKGROUND:
            print('subtracting mean')
            usedstack = stack[frames, :, :] - stack[:, :, :].mean(axis=0) # remove space background

        elif method == DIFF_WITH_ZERO:
            print('subtracting frame 0')
            usedstack = stack[frames, :, :] - stack[0, :, :]

        elif method == DIFF_WITH_PREVIOUS:
            print('subtracting previous frame')
            usedstack = stack[np.array(frames[:-1])+1, :, :] - stack[frames[:-1], :, :]

        elif method == NONE:
            usedstack = stack[frames, :, :]

        elif method == TWOCHANNEL:
            usedstack = np.zeros((len(frames), stacks[0].shape[1], stacks[0].shape[2], 3))
            for i, this_stack in enumerate(stacks):
                print(this_stack.shape)
                if stackcolors[i] == 'red':
                    usedstack[:, :, :, 0] = normalise_stack(this_stack[frames, :, :], file+' red')

                elif stackcolors[i] == 'green':
                    usedstack[:, :, :, 1] = normalise_stack(this_stack[frames, :, :], file+' green')
                else:
                    assert False, 'colors other than red and green not yet implemented'


    # common.term_hist(usedstack)

    # print('min mean max', usedstack.min(), usedstack.mean(), usedstack.max())
    if not no_stack:
        vmin = np.quantile(usedstack, 0.01)
        vmax = np.quantile(usedstack, 0.99)

    def show(index):
        if backwards:
            timestep = frames[-index]
        else:
            timestep = frames[index]

        ax.clear()

        if not no_stack:

            if len(usedstack.shape) == 4:
                frame = usedstack[index, :, :, :]
            else:
                frame = usedstack[index, :, :]

            color = 'white' if frame.mean()/(frame.max()-frame.min()) < 0.2 else 'black' # this used to be usedstack not frame
        
            show_single_frame(ax, frame, pixel_size, channel=channel, vmin=vmin, vmax=vmax, window_size_x=window_size_x, window_size_y=window_size_y)

        else:
            color = 'gray'

            show_single_frame(ax, None, pixel_size, channel=channel, window_size_x=window_size_x, window_size_y=window_size_y)
        
        time_string = speed_string(time_mult, every_nth_frame*nth_frame)#+f'\ntime = {timestep*time_step*nth_frame:.1f}s'
        ax.text(0.95, 0.05, time_string, color=color, transform=ax.transAxes, ha='right', fontsize=10)

        ax.text(0.1, 0.9, dataset_name, transform=ax.transAxes, fontsize=15, color=color)
     
        
        # for hiding border
        fig.set_size_inches(*fig.get_size_inches()) # apparently this is needed to make subplots_adjust work
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

    
        func(timestep, ax)

    if no_stack:
        num_timesteps = max_num_frames
        assert num_timesteps > 0
    else:
        num_timesteps = usedstack.shape[0]

    common.save_gif(show, range(num_timesteps), fig, outputfilename, fps=fps, dpi=dpi)

def show_single_frame(ax, frame, pixel_size, window_size_x, window_size_y, channel=None, vmin=None, vmax=None):

    
    # vmin = np.quantile(usedstack, 0.01)
    # vmax = np.quantile(usedstack, 0.99)
    # vmin = stack.min()
    # vmax = stack.max()
    # vmin = -5.5e-2
    # vmax = 6.2e-2
    # vmin = None
    # vmax = None

    if channel == 'red':
        colors = [(0, 0, 0), (1, 0, 0)] # black > red
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("Custom", colors, N=20)
    elif channel == 'green':
        colors = [(0, 0, 0), (0, 1, 0)] # black > green
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("Custom", colors, N=20)
    else:
        cmap = matplotlib.cm.Greys

    if frame != None:
        ax.imshow(frame, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, extent=(0, window_size_x, 0, window_size_y))
    
    ax.set_xlim(0, window_size_x)
    ax.set_ylim(0, window_size_y)
    
    if frame:
        color = 'white' if frame.mean()/(frame.max()-frame.min()) < 0.2 else 'black'
    else:
        color = 'gray'
    common.add_scale_bar(ax, pixel_size, color=color)
   


def go(file, outputfilename, add_drift=False, display_small=False, method=NONE, highlights=False, backward=False, flip_y=False):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        stack      = data['stack']
        pixel_size = data['pixel_size']
        time_step  = data['time_step']
        channel    = data.get('channel')
        window_size_x = data['window_size_x']
        window_size_y = data['window_size_y']

        if add_drift:
            stack = common.add_drift_intensity(stack, 1)
            print('adding drift')
            stack = np.swapaxes(stack, 1, 2)

        if display_small:
            # crop
            stack = stack[:, :300, :300]

        if flip_y:
            stack = stack[:, ::-1, :]

        # print(stack.shape[1], 'x', stack.shape[2], 'px')
        # print(stack.min(), stack.max())


        save_array_movie(stack, pixel_size, time_step, file, outputfilename,
                         nth_frame=data.get('nth_frame', 1),
                         display_small=DISPLAY_SMALL, inverse_colors=INVERSE_COLORS, highlights=HIGHLIGHTS,
                         backwards=BACKWARDS, method=METHOD, channel=channel,
                         dataset_name=data.get('NAME'),
                         window_size_x=window_size_x, window_size_y=window_size_y
        )
        # save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cmd31/figures/{filename}.mp4",
        #                  nth_frame=data.get('nth_frame', 1))
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4")

        
if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        
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
        if ADD_DRIFT:
            filename += '_drifted'

        go(
            file,
            outputfilename=filename,
            add_drift=ADD_DRIFT,
            method=METHOD,
            display_small=DISPLAY_SMALL,
            highlights=HIGHLIGHTS,
            backward=BACKWARDS
        )