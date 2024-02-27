import numpy as np
import common

stack = common.load('raw_data/0.34_s2_s3.npy')
np.savez('preprocessing/data/stack_eleanor.npz', stack=stack, time_step=0.5, pixel_size=0.288, particle_diameter=2.82)