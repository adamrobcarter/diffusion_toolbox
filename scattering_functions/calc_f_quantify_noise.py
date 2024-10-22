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

######## first
# scattering_functions.calc_both.calc_fnor_f_type(
#     'F',
#     log=False,
#     # max_K=1.08,
#     num_k_bins=50,  # computation time is proportional to this squared
#     # file_suffix='_25bins',
#     cores=16, # increasing this above 16 seems risky
#     max_time_origins=10000, # all time origins
#     first_point_only=True,
#     file_suffix='_first',
# )

####### for testing k_x=0
scattering_functions.calc_both.calc_for_f_type(
    'F',
    log=False,
    # max_K=1.08,
    num_k_bins=50,  # computation time is proportional to this squared
    # file_suffix='_25bins',
    cores=16, # increasing this above 16 seems risky
    max_time_origins=50, # computation time is directly proportional # eleanorlong needs this big
    # S_only=True,
    use_zero=True,
    file_suffix='_usezero',
)