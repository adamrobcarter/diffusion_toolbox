import DDM.phiDM
import common
import numpy as np
import time

background_removed = True

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']
    time_step  = data['time_step']

    t0 = time.time()
    DDM.phiDM.calc(stack, pixel_size)
    t1 = time.time()

    # np.savez(f'DDM/data/ddm_{file}.npz', k=k, F_D_sq=F_D_sq, t=t,
    #          background_removed=background_removed, use_every_nth_frame=use_every_nth_frame,
    #          computation_time_ddm=t1-t0, num_k_bins=num_k_bins,
    #          pixel_size=pixel_size, particle_diameter=data['particle_diameter'])