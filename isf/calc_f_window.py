import isf.calc_both
import common
"""

OUT OF DATE, SHOULD USE window PARAM OF calc_f.py

# doubleside
for file in common.files_from_argv('particle_detection/data', 'particles_'):
    d_frames = None

    if '001' in file:
        max_time_origins = 10000
    elif '010' in file:
        max_time_origins = 1000
    elif '034' in file:
        max_time_origins = 300
    else:
        max_time_origins = 600

    max_time_origins /= 2
    max_time_origins /= 2
    # if 'sim_nohydro' in file:
    #     max_time_origins /= 4
    # if 'cropsquare' in file:
    #     max_time_origins *= 2

    if file.startswith('sim_nohydro'):
        max_time_origins = 250

    if file == 'sim_nohydro_011_L320_trim7.11hr':
        max_time_origins = 4000
        d_frames = [0, 1, 16]

    if file == 'sim_nohydro_011_L640':
        max_time_origins = 4000
        # d_frames = [0, 64,]
        d_frames = [0, 1, 4, 16, 64,]

    d_frames = [0, 1, 4, 16, 64, 256, 1024, 4096]

    if '_pot_longer' in file:
        d_frames = [0, 1, 4, 16, 64, 256, 1024]

    window = True
    
    isf.calc_both.calc_for_f_type(
        file,
        'F',
        cores=16, # increasing this above 16 seems risky for big datasets
        max_time_origins=max_time_origins, # computation time is directly proportional
        use_zero=True, use_doublesided_k=False,
        # file_suffix='',
        d_frames=d_frames,
        num_k_bins=60,
        window = 'blackmanharris' if window else None,
        file_suffix = '_bhwindow' if window else ''
    )

    window = False
    
    isf.calc_both.calc_for_f_type(
        file,
        'F',
        cores=16, # increasing this above 16 seems risky for big datasets
        max_time_origins=max_time_origins, # computation time is directly proportional
        use_zero=True, use_doublesided_k=False,
        # file_suffix='',
        d_frames=d_frames,
        num_k_bins=60,
        window = 'blackmanharris' if window else None,
        file_suffix = '_bhwindow' if window else ''
    )