import DDM.DDM
import DDM.show
import common
import time
import DDM.viz

background_removed = False

for file in common.files_from_argv('preprocessing/data', 'stack_'):
    data = common.load(f'preprocessing/data/stack_{file}.npz')
    stack      = data['stack']
    pixel_size = data['pixel_size']
    time_step  = data['time_step']
    NAME       = data.get('NAME')
    channel    = data.get('channel')

    if background_removed:
        stack = stack - stack.mean(axis=0)

    num_k_bins = 100

    # function to save every 10 iterations of the calculation and also at the beginning since we are gready
    def callback(i, k, F_D_sq, F_D_sq_unc, t, F_D_sq_all, time_origins, F_D_sq_noradial, k_x, k_y):
        if i % 10 == 0 or i in [1, 3, 5, 15]:
            sigma = 0
            pixel = 0
            
            # DDM.show.show(file, k, F_D_sq, F_D_sq_unc, t, sigma, pixel, live=True, NAME=NAME, channel=channel)
            # DDM.viz.show(file, F_D_sq_all, t, k, time_origins)

            common.save_data(f'DDM/data/ddm_{file}.npz', quiet=True,
                k=k, F_D_sq=F_D_sq, F_D_sq_unc=F_D_sq_unc, t=t,
                background_removed=background_removed, use_every_nth_frame=0,
                num_k_bins=num_k_bins,
                F_D_sq_all=F_D_sq_all, time_origins=time_origins, F_D_sq_noradial=F_D_sq_noradial, k_x=k_x, k_y=k_y,
                pixel_size=pixel_size, particle_diameter=data.get('particle_diameter'), NAME=data.get('NAME'),
                pack_frac_given=data.get('pack_frac_given'), partial=i)

    # do the DDM calculation
    t0 = time.time()
    k, t, F_D_sq, F_D_sq_unc, use_every_nth_frame, F_D_sq_all, time_origins,F_D_sq_noradial, k_x, k_y = DDM.DDM.calc(stack, pixel_size, time_step, num_k_bins, callback=callback)
    t1 = time.time()

    common.save_data(f'DDM/data/ddm_{file}.npz', k=k, F_D_sq=F_D_sq, F_D_sq_unc=F_D_sq_unc, t=t,
        background_removed=background_removed, use_every_nth_frame=use_every_nth_frame,
        computation_time_ddm=t1-t0, num_k_bins=num_k_bins,
        F_D_sq_all=F_D_sq_all, time_origins=time_origins, F_D_sq_noradial=F_D_sq_noradial, k_x=k_x, k_y=k_y,
        pixel_size=pixel_size, particle_diameter=data.get('particle_diameter'), NAME=data.get('NAME'),
        pack_frac_given=data.get('pack_frac_given'), complete=True)