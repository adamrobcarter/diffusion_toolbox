import isf.calc_both
import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):

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
            max_time_origins *= 16

        if file == 'sim_nohydro_010_L544' or file == 'brennan_hydro_010_L544':
            max_time_origins *= 4

        cores = 16 # increasing this above 16 seems risky. if the program freezes, try reducing this
        if 'L1280' in file:
            cores = 4

        for max_time_origins in [250, 1000, 4000, 16000]:
            isf.calc_both.calc_for_f_type(
                file,
                'F',
                num_k_bins=30,
                linear_log_crossover_k=0.5,
                # file_suffix='_25bins',
                cores=cores,
                max_time_origins=max_time_origins,
                d_frames=[0, 1],
                file_prefix = 'first_',
                file_suffix = f'_mto{max_time_origins}'
            )
