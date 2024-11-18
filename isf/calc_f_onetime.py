import isf.calc_both
import common

# doublesided
for file in common.files_from_argv('particle_detection/data', 'particles_'):
    if '001' in file:
        max_time_origins = 10000
    elif '010' in file:
        max_time_origins = 1000
    elif '034' in file:
        max_time_origins = 300
    else:
        max_time_origins = 150

    # max_time_origins /= 2
    # if 'cropsquare' in file:
    #     max_time_origins *= 2
    
    isf.calc_both.calc_for_f_type(
        file,
        'F',
        cores=16, # increasing this above 16 seems risky for big datasets
        max_time_origins=max_time_origins, # computation time is directly proportional # eleanorlong needs this big (~1000)
        use_zero=True, use_doublesided_k=True,
        file_suffix='_doublesided_onetime',
        d_frames=[0, 120],
        num_k_bins=20,
    )