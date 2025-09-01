import isf.calc_both
import common

def go(file, quiet=False):
    if '001' in file:
        max_time_origins = 10000
    elif '002' in file:
        max_time_origins = 5000
    elif '010' in file:
        max_time_origins = 1000
    elif '034' in file:
        max_time_origins = 300
    else:
        max_time_origins = 150
    
    # no more changing! if you need to change for a specific file, add an if

    max_time_origins *= 2

    if file.startswith('sim_nointer'):
        # max_time_origins *= 16
        pass
        
    if 'hydro_011_L' in file or 'brennan_hydro_011_L' in file:
        max_time_origins = 4000

    cores = 16 # increasing this above 16 seems risky. if the program freezes, try reducing this
    if 'L1280' in file:
        cores = 4

    # if file == 'brennan_hydro_034':

    if 'div64' in file:
        max_time_origins /= 64
        

    # print('DEBUG')
    # max_time_origins /= 1000

    # if 'eleanor' in file:
    #     window = 'blackmanharris'
    #     max_time_origins *= 1 / 0.132
    # else:
    #     window = None

    d_frames = [0, 0.5] if 'mixt' in file else [0, 1]

    if file == 'sim_nohydro_011_L320_test_singlet': d_frames = [0, 16]

    isf.calc_both.calc_for_f_type(
        file,
        'F',
        num_k_bins=60,
        # file_suffix='_25bins',
        cores=cores,
        max_time_origins=max_time_origins,
        d_frames=d_frames,
        file_prefix = 'first_',
        # file_prefix = 'first32_'
        # window=window,
        # file_suffix='_new',
        quiet=quiet
    )

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        go(file)