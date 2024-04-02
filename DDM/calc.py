import DDM.DDM
import common
import numpy as np

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']
    time_step  = data['time_step']

    k, t, F_D = DDM.DDM.calc(stack, pixel_size, time_step)
    np.savez(f'DDM/data/ddm_{file}.npz', k=k, F_D=F_D, t=t)