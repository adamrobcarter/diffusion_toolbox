import common
import numpy as np
import time
import tqdm

WINDOW_LENGTH = 5

t0 = time.time()

if __name__ == '__main__':
    for file in common.files_from_argv('preprocessing/data/', 'stack_'):
        data = common.load(f'preprocessing/data/stack_{file}.npz')
        stack = data['stack']
        width = stack.shape[1]
        height = stack.shape[2]

        # we do it in place

        for i in tqdm.trange(stack.shape[0]-WINDOW_LENGTH):
            stack[i] = stack[i:i+WINDOW_LENGTH, :, :].mean(axis=0)

        stack = stack[:-WINDOW_LENGTH, :, :] # remove excess


        common.save_data(f'preprocessing/data/stack_{file}_rolling{WINDOW_LENGTH}.npz',
                        stack=stack, time_step=data['time_step'], pixel_size=data['pixel_size'],
                        pack_frac_given=data.get('pack_frac_given'), particle_diameter=data.get('particle_diameter'),
                        nth_frame=data.get('nth_frame', 1))
        
        nth_frame = (stack.shape[0]+WINDOW_LENGTH)//50 # the +WINDOW_LENGTH is so that we use the same frames for nth_frame is if we didn't do the rolling av, to allow comparisons
        common.save_data(f'preprocessing/data/stack_{file}_rolling{WINDOW_LENGTH}_small.npz',
                        stack=stack[::nth_frame, :, :], time_step=data['time_step']*nth_frame, pixel_size=data['pixel_size'],
                        pack_frac_given=data.get('pack_frac_given'), particle_diameter=data.get('particle_diameter'),
                        nth_frame=data.get('nth_frame', 1)*nth_frame)
        

