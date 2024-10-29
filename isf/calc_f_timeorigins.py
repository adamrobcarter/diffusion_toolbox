import scattering_functions.calc_both

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

# for max_time_origins in [1, 2, 4, 8, 16, 32, 64]:
for max_time_origins in [64]:

    scattering_functions.calc_both.calc_for_f_type(
        'F',
        log=False,
        # max_K=1.08,
        num_k_bins=20,  # computation time is proportional to this squared
        file_suffix=f'_timeorigins{max_time_origins}',
        cores=16, # increasing this above 16 seems risky (16 times the RAM)
        max_time_origins=max_time_origins, # computation time is directly proportional
        # S_only=True,
    )