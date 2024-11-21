import isf.calc_both
import common

# scattering_functions.calc_both.calc_for_f_type('F')
# scattering_functions.calc_both.calc_for_f_type(
#     'F',
#     log=False,
#     max_K=0.35,
#     num_k_bins=200,  # computation time is proportional to this squared
#     file_suffix='_200bins',
#     cores=60, 
#     max_time_origins=10, # computation time is directly proportional
#     S_only=True
# )

# scattering_functions.calc_both.calc_for_f_type(
#     'F',
#     log=False,
#     max_K=0.35,
#     num_k_bins=100,  # computation time is proportional to this squared
#     file_suffix='_100bins',
#     cores=60, 
#     max_time_origins=10, # computation time is directly proportional
#     S_only=True
# )

# scattering_functions.calc_both.calc_for_f_type(
#     'F',
#     log=False,
#     max_K=0.35,
#     num_k_bins=50,  # computation time is proportional to this squared
#     file_suffix='_50bins',
#     cores=60, 
#     max_time_origins=10, # computation time is directly proportional
#     S_only=True
# )

####### normal
# scattering_functions.calc_both.calc_for_f_type(
#     'F',
#     log=False,
#     # max_K=1.08,
#     num_k_bins=50,  # computation time is proportional to this squared
#     # file_suffix='_25bins',
#     cores=16, # increasing this above 16 seems risky
#     max_time_origins=500, # computation time is directly proportional # eleanorlong needs this big
#     # S_only=True,
# )

####### for testing k_x=0
# for file in common.files_from_argv('particle_detection/data', 'particles_'):
#     scattering_functions.calc_both.calc_for_f_type(
#         file,
#         'F',
#         log=False,
#         # max_K=1.08,
#         num_k_bins=50,  # computation time is proportional to this squared
#         # file_suffix='_25bins',
#         cores=16, # increasing this above 16 seems risky
#         max_time_origins=25, # computation time is directly proportional # eleanorlong needs this big
#         # S_only=True,
#         use_zero=True,
#         file_suffix='_usezero',
#     )


# for file in common.files_from_argv('particle_detection/data', 'particles_'):
#     scattering_functions.calc_both.calc_for_f_type(
#         file,
#         'F',
#         log=False,
#         # max_K=1.08,
#         # file_suffix='_25bins',
#         cores=16, # increasing this above 16 seems risky
#         max_time_origins=10000, # computation time is directly proportional # eleanorlong needs this big
#         # S_only=True,
#         use_zero=False,
#         # file_suffix='_smallbins',
#     )

# for eleanorlong034_crop
# for file in common.files_from_argv('particle_detection/data', 'particles_'):
#     scattering_functions.calc_both.calc_for_f_type(
#         file,
#         'F',
#         log=False,
#         # max_K=1.08,
#         # file_suffix='_25bins',
#         cores=16, # increasing this above 16 seems risky
#         max_time_origins=10000, # computation time is directly proportional # eleanorlong needs this big
#         # 10000 is very smooth for el034
#         # el001 try 50000
#         # S_only=True,
#         use_zero=False,
#         file_suffix='_first',
#         d_frames=[0, 1],
#         num_k_bins=20
#     )
# you need to call the nocrop one crop1 in order to stop it getting modified by something else

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

    max_time_origins /= 2
    max_time_origins /= 2
    # if 'sim_nohydro' in file:
    #     max_time_origins /= 4
    # if 'cropsquare' in file:
    #     max_time_origins *= 2

    if file.startswith('sim_nohydro_011_L320'):
        max_time_origins = 1000
    
    isf.calc_both.calc_for_f_type(
        file,
        'F',
        cores=16, # increasing this above 16 seems risky for big datasets
        max_time_origins=max_time_origins, # computation time is directly proportional
        use_zero=True, use_doublesided_k=False,
        # file_suffix='',
        # d_frames=[0, 120],
        num_k_bins=30,
        linear_log_crossover_k=0.5
    )