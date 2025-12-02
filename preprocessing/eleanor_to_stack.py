import numpy as np
import common

if __name__ == '__main__':
    # stack = common.load('raw_data/0.34_s2_s3.npy')
    # common.save_data('preprocessing/data/stack_eleanor0.34.npz', stack=stack, time_step=0.5, 
    #          pixel_size=0.288, particle_diameter=2.82, pack_frac_given=0.34)

    stack = common.load('raw_data/0.01_s2_s6.npy')
    stack = stack.transpose([2, 0, 1])
    # stack = stack[:stack.shape[0]//4, :, :] WHY WAS THIS HERE?
    common.save_data('preprocessing/data/stack_eleanor0.01.npz', stack=stack, time_step=0.5, 
            pixel_size=0.288, particle_diameter=2.82, pack_frac_given=0.01)
    common.save_data('preprocessing/data/stack_eleanor0.01_first.npz', stack=stack[[0], :, :], time_step=0.5, 
            pixel_size=0.288, particle_diameter=2.82, pack_frac_given=0.01)