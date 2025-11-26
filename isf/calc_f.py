import isf.calc_both
import common

def go(file, window=False):
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

    d_frames = [0, 1, 4, 16, 64, 256, 1024]
    # d_frames = [0, 1, 4, 16, 64, 256, 1024]

    if '_pot_longer' in file or file=='eleanor0.34':
        d_frames = [0, 1, 4, 16, 64, 256, 1024]

    # if 'L2560' in file:
    if '_32' in file:
        d_frames = [0, 1, 2, 8, 32, 128, 512]
    else:
        d_frames = [0, 1, 4, 16, 64, 256]
    # if 'nbody_open' in file:
    #     d_frames = [0, 1, 4, 16, 64]

    if '_1s' in file:
        d_frames = [0, 1]

    if file == 'ld_hydro_dpstokes_0.114_L1280_short':
        d_frames = [0, 2]

    # d_frames = common.exponential_integers(0, 4096)
    # d_frames = d_frames[d_frames < 623]
    # print("NOVTE 60/2 below")
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

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file)