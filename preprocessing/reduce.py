import common
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize

FIRST_FRAME = False
FRAME_DIFF = True
REMOVE_BKG = False

def do_reduce():

    for file in common.files_from_argv('preprocessing/data', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack      = data['stack']
        pixel_size = data['pixel_size']
        time_step = data['time_step']

    while True:

        print(f'Total size is {stack.shape[0]} frames in time')
        first_frame = int(input(f'Please enter FIRST frame you want (btwn 1 and {stack.shape[0]}): '))
        last_frame = int(input(f'Please enter LAST frame you want (btwn 1 and {stack.shape[0]}): '))
        interval = int(input(f'Please enter every nth frame to use: '))
        
        print(f'Will make a reduced stack from {first_frame} to {last_frame}')
        stackr = stack[first_frame-1:last_frame-1:interval, :, :]
        common.save_data(f'preprocessing/data/stack_{file}_red_{first_frame}_{last_frame}_{interval}.npz', stack=stackr, pixel_size=pixel_size, time_step=time_step*interval, nth_frame=interval)
        del stackr
        print(f'Saved a reduced stack as {file}_red_{first_frame}_{last_frame}_{interval}, Call this in the future')
    
        inputted = input(f'Do you want to do another stack? if YES enter 1 otherwise return')
        if not inputted:
            break

if __name__ == '__main__':

    do_reduce()