import numpy as np
import matplotlib.pyplot as plt
import common
import sys
import matplotlib.cm
import warnings

DISPLAY_SMALL = False
INVERSE_COLORS = False
SHOW_TIMESTEP = True
HIDE_ANNOTATIONS = False

SHOW_HISTOGRAMS = False

HIGHLIGHTS = True # displays HIGHLIGHTS_NUMBER frames evenly throughout the stack instead of the first 50
HIGHLIGHTS_NUMBER = 10
# HIGHLIGHTS = True
# BACKWARDS = True
BACKWARDS = False # plays the stack backwards. DIFF_WITH_ZERO is now compared to the last frame

# possible display methods
REMOVE_BACKGROUND = 1
DIFF_WITH_ZERO = 2 # frame - frame1
DIFF_WITH_PREVIOUS = 3
NONE = 4
TWOCHANNEL = 5
DIFF_WITH_TEN = 6
FIRSTLAST = 7

METHOD_NAME = {
    REMOVE_BACKGROUND: 'i(t) - mean(i)',
    DIFF_WITH_ZERO: 'i(t) - i(0)',
    DIFF_WITH_PREVIOUS: 'i(t) - i(t-1)',
    DIFF_WITH_TEN: 'i(t) - i(t-10)',
    NONE: '',
    TWOCHANNEL: '',
}

# select your method here
# METHOD = DIFF_WITH_PREVIOUS
# METHOD = DIFF_WITH_ZERO
# METHODS = [NONE, DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_PREVIOUS]
# METHODS = [NONE, DIFF_WITH_TEN,REMOVE_BACKGROUND, DIFF_WITH_ZERO, DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_TEN]
METHODS = [NONE]
# METHODS = [NONE, DIFF_WITH_PREVIOUS, REMOVE_BACKGROUND]
# METHODS = [DIFF_WITH_PREVIOUS]
# METHODS = [DIFF_WITH_ZERO]
# METHODS = [REMOVE_BACKGROUND]
# METHODS = [FIRSTLAST]

ADD_DRIFT = False

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
                     func=lambda timestep, ax : None, nth_frame=1, max_num_frames=30,
                     dpi=300,
                     display_small=True, inverse_colors=False, highlights=False,
                     backwards=False, method=NONE, stacks=None, stackcolors=None, channel=None,
                     dataset_name=None,
                     num_timesteps_in_data=None, # how many timesteps are in the original data?. only needed if stack is None
                     output_type='movie',
                     figsize_mult=1.0,
                     every_nth_frame=None,
                     annotation_color=None,
                     ):
    if method == NONE:
        assert not backwards


    no_stack = stack is None

    if no_stack:
        assert num_timesteps_in_data > 0, f'num_timesteps_in_data = {num_timesteps_in_data}'
        figsize = np.array([3, 3])
    else:
        if method == TWOCHANNEL:
            assert stack is None
            assert stacks is not None
            stack = stacks[0] # this copy is so that we can use the same checks in twochannel mode as we do in one channnel mode
        else:
            assert stack is not None
            assert stacks is None

        num_timesteps_in_data = stack.shape[0]
        assert num_timesteps_in_data > 0, f'num_timesteps_in_data = {num_timesteps_in_data}'
    
        figsize = np.array(stack.shape)[[2, 1]] / dpi
        while figsize.mean() < 3:
            figsize *= 2
        # if display_small and not : # potentially remove this for marine
        #     figsize = (2.3, 2.3)

    figsize = figsize * figsize_mult
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # stack = np.transpose(stack, [0, 2, 1]) # idk why, check on clearly rectangular data like psiche089_slice0


    time_mult = 1
    if file.startswith('marine'):
        time_mult = 0.25
    if file.startswith('marine'):
        time_mult = 0.25
    
    if every_nth_frame is None:
        if highlights:
            print('HIGHLIGTHS')
            every_nth_frame = num_timesteps_in_data // HIGHLIGHTS_NUMBER
            max_num_frames = HIGHLIGHTS_NUMBER
        else:
            every_nth_frame = 1

    if file == 'pierre_exp':
        every_nth_frame = 10
        time_mult = 0.1

    assert every_nth_frame > 0, f'every_nth_frame = {every_nth_frame}, is your data too short?'

    fps = 1/time_step * time_mult

    if method == FIRSTLAST:
        fps = 1

    # these are cause of the reduce operation
    if file.endswith('_25'):
        fps *= 25
    fps *= nth_frame#*every_nth_frame


    if fps < 0.2:
        warnings.warn(f'fps = {fps}')

    if not no_stack:
        num_timesteps_in_data = stack.shape[0]
        assert num_timesteps_in_data > 0, f'num_timesteps_in_data = {num_timesteps_in_data}'
    
    frames = range(0, min(num_timesteps_in_data, max_num_frames*every_nth_frame), every_nth_frame)

    assert len(frames) > 0, f'frames = {frames}, num_timesteps_in_data = {num_timesteps_in_data}, max_num_frames = {max_num_frames}, every_nth_frame = {every_nth_frame}'

    if (not no_stack) and SHOW_HISTOGRAMS:
        print('stack[frames]:')
        common.term_hist(stack[frames])

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

        elif method == DIFF_WITH_TEN:
            print('subtracting previous frame')
            usedstack = stack[np.array(frames[:-1])+10, :, :] - stack[np.array(frames[:-1])+0, :, :]
            print('I want to look later in the stack')

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

        elif method == FIRSTLAST:
            print(f'delta t = {time_step * stack.shape[0]}s')
            usedstack = stack[[0, -1], :, :]

    if (not no_stack) and SHOW_HISTOGRAMS:
        print('used_stack:')
        common.term_hist(usedstack)

    # print('min mean max', usedstack.min(), usedstack.mean(), usedstack.max())
    if not no_stack:
        vmin = np.quantile(usedstack, 0.01)
        vmax = np.quantile(usedstack, 0.99)

    print('frames', frames)

    def show(index):
        if backwards:
            timestep = frames[-index]
        else:
            timestep = frames[index]

        ax.clear()

        if method == FIRSTLAST:
            # start, end = ax.get_xlim()
            print('set_ticks')
            # ax.xaxis.set_ticks(np.arange(start, end, 10))
            import matplotlib.ticker as plticker
            loc = plticker.MultipleLocator(base=10) # this locator puts ticks at regular intervals
            ax.xaxis.set_major_locator(loc)


        if not no_stack:

            if len(usedstack.shape) == 4:
                frame = usedstack[index, :, :, :]
            else:
                frame = usedstack[index, :, :]

            # print(frame.mean()/(frame.max()-frame.min()))
            if annotation_color:
                color = annotation_color
            else:
                color = 'white' if frame.mean()/(frame.max()-frame.min()) > 0.2 else 'black' # this used to be usedstack not frame
                                                                                         # also the > used to be a <,
        
            show_single_frame(file, ax, frame, pixel_size, channel=channel, vmin=vmin, vmax=vmax,
                              window_size_x=window_size_x, window_size_y=window_size_y,
                              hide_scale_bar=HIDE_ANNOTATIONS, annotation_color=color)

        else:
            color = 'gray'

            show_single_frame(file, ax, None, pixel_size, channel=channel, window_size_x=window_size_x, window_size_y=window_size_y,
                              hide_scale_bar=HIDE_ANNOTATIONS)
        
        time_string = speed_string(time_mult, every_nth_frame*nth_frame)
        if not HIDE_ANNOTATIONS:
            if SHOW_TIMESTEP:
                time_string = ''
                # time_string += f'\nframe = {timestep*nth_frame:.0f}'
                time_string += f'\nt = {timestep*time_step*nth_frame:.1f}s'
                ax.text(0.95, 0.05, time_string, color=color, transform=ax.transAxes, ha='right', fontsize=10)

            ax.text(0.05, 0.9, dataset_name, transform=ax.transAxes, fontsize=15, color=color)
            
            # ax.text(0.05, 0.82, METHOD_NAME[METHOD], transform=ax.transAxes, fontsize=15, color=color)
     
        
        # for hiding border
        if output_type == 'frames':
            # this used to be commented out, it seems we don't need it for save_gif
            # but we might need it for save_frames
            fig.set_size_inches(*fig.get_size_inches()) # apparently this is needed to make subplots_adjust work
            fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

    
        func(timestep, ax)

    if no_stack:
        # num_timesteps = max_num_frames # cba to work out why this didn't work
        # assert num_timesteps > 0
        num_timesteps = len(frames)
        assert num_timesteps > 0, f'num_timesteps = len(frames) = {num_timesteps}'
    else:
        num_timesteps = usedstack.shape[0]
        assert num_timesteps > 0, f'num_timesteps = usedstack.shape[0] = {num_timesteps}'

    if output_type == 'movie':
        common.save_gif(show, range(num_timesteps), fig, outputfilename, fps=fps, dpi=dpi)
    elif output_type == 'frames':
        common.save_frames(show, range(num_timesteps), fig, outputfilename, file, fps, dpi=dpi)
    else:
        raise ValueError(f'output must be "movie" or "frames", not {output_type}')

def show_single_frame(file, ax, frame, pixel_size, window_size_x, window_size_y, channel=None, vmin=None, vmax=None,
                      hide_scale_bar=False, annotation_color='black'):

    
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

    if frame is not None:
        if 'faxtor' in file:
            cmap = matplotlib.cm.Greys_r
        # ax.imshow(frame, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, extent=(0, window_size_x, 0, window_size_y))
        X, Y = np.meshgrid(np.linspace(0, frame.shape[0], frame.shape[0])*pixel_size, np.linspace(0, frame.shape[1], frame.shape[1])*pixel_size, indexing='ij')
        # X, Y = np.meshgrid(np.arange(0, window_size_x, pixel_size), np.arange(0, window_size_y, pixel_size), indexing='ij')
        ax.pcolormesh(X, Y, frame, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')

    if 'unwrap' in file:
        ax.set_xlim(-window_size_x, 2*window_size_x)
        ax.set_ylim(-window_size_y, 2*window_size_y)
    else:
        ax.set_xlim(0, window_size_x)
        ax.set_ylim(0, window_size_y)
    # ax.set_xlim(-window_size_x, window_size_x)
    # ax.set_ylim(-window_size_x, window_size_y)
    ax.set_aspect('equal')

    if not hide_scale_bar:
        common.add_scale_bar(ax, pixel_size, color=annotation_color)
   


def go(file, data=None, outputfilename='please provide me', add_drift=False, display_small=False, method=NONE, highlights=False, backward=False, flip_y=False,
       output_type='movie', figsize_mult=1.0, dpi=200, crop=None):
        
        if data is None:
            # we allow the user to pass in data for legacy reasons 
            # (it was so they could generate movies with different methods without reloading the data from disk)
            data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        stack      = data['stack']
        pixel_size = data['pixel_size']
        time_step  = data['time_step']
        channel    = data.get('channel')

        print('time:')
        if SHOW_HISTOGRAMS:
            common.term_bar(stack[:, :, :].mean(axis=(1, 2)), range(0, stack.shape[0]+1))

        if add_drift:
            stack = common.add_drift_intensity(stack, 1)
            print('adding drift')
            # stack = np.swapaxes(stack, 1, 2) why tf was this here?

        if crop:
            stack = stack[:, crop[0][0]:crop[0][1], crop[1][0]:crop[1][1]]
            assert stack.size

        if display_small:
            crop = 250
            x_start = int((stack.shape[1] - crop)/2)
            y_start = int((stack.shape[2] - crop)/2)
            stack = stack[:, x_start:x_start+crop, y_start:y_start+crop]

        if flip_y:
            stack = stack[:, ::-1, :]
            
        window_size_x = stack.shape[1] * pixel_size
        window_size_y = stack.shape[2] * pixel_size

        # print(stack.shape[1], 'x', stack.shape[2], 'px')
        print('stack min nmax', stack.min(), stack.max())


        save_array_movie(stack, pixel_size, time_step, file, outputfilename,
                         nth_frame=data.get('nth_frame', 1),
                         display_small=DISPLAY_SMALL, inverse_colors=INVERSE_COLORS, highlights=highlights,
                         backwards=backward, method=method, channel=channel,
                         dataset_name=data.get('NAME'),
                         window_size_x=window_size_x, window_size_y=window_size_y,
                         output_type=output_type,
                         figsize_mult=figsize_mult,
                         dpi = dpi,
        )
        # save_array_movie(stack, pixel_size, time_step, file, f"/home/acarter/presentations/cmd31/figures/{filename}.mp4",
        #                  nth_frame=data.get('nth_frame', 1))
        # save_array_movie(stack_copy, pixel_size, time_step, file, f"/home/acarter/presentations/cin_first/figures/{filename}.mp4")

        
if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        # num = int(file[6:9])
        # if num < 114:
        #     continue
        
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        
        # try:
        if True:
            for METHOD in METHODS:
                for crop in [False]:
                # for crop in [False, True]:
            
                    filename = f'preprocessing/figures_png/stack_movie_{file}'

                    if crop:
                        filename += '_crop'

                    if METHOD == REMOVE_BACKGROUND:
                        filename += '_bkgrem'
                    elif METHOD == DIFF_WITH_ZERO:
                        filename += '_diffframe1'
                    elif METHOD == DIFF_WITH_PREVIOUS:
                        filename += '_diffprev'
                    elif METHOD == DIFF_WITH_TEN:
                        filename += '_diff10'
                    if METHOD == FIRSTLAST:
                        filename += '_firstlast'
                    if HIGHLIGHTS:
                        filename += '_highlights'
                    if BACKWARDS:
                        filename += '_backwards'
                    if ADD_DRIFT:
                        filename += '_drifted'
                    
                    filename += '.gif'

                    go(
                        file,
                        data,
                        outputfilename=filename,
                        add_drift=ADD_DRIFT,
                        method=METHOD,
                        display_small=crop,
                        highlights=HIGHLIGHTS,
                        backward=BACKWARDS,
                        output_type='movie',
                    )
        # except Exception as err:
        #     print()
        #     print(err)
        #     print()