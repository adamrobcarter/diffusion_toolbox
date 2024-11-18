import DDM.DDM as DDM
import common
import numpy as np
import time

background_removed = False

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']
    time_step  = data['time_step']

    if background_removed:
        stack = stack - stack.mean(axis=0)

    t0 = time.time()
    num_k_bins = 50
    k, t, F_D_sq, F_D_sq_unc, use_every_nth_frame, F_D_sq_all, time_origins = DDM.calc(stack, pixel_size, time_step, num_k_bins)
    t1 = time.time()

    common.save_data(f'DDM/data/ddm_{file}.npz', k=k, F_D_sq=F_D_sq, F_D_sq_unc=F_D_sq_unc, t=t,
        background_removed=background_removed, use_every_nth_frame=use_every_nth_frame,
        computation_time_ddm=t1-t0, num_k_bins=num_k_bins,
        pixel_size=pixel_size, particle_diameter=data.get('particle_diameter'),
        pack_frac_given=data.get('pack_frac_given'),
        NAME=data.get('NAME'), channel=data.get('channel'))