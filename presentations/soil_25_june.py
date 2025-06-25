import preprocessing.stack_movie
import matplotlib.pyplot as plt
import DDM.show_single
import common

path = '/home/acarter/presentations/soil_25_june/figures'

# preprocessing.stack_movie.go(
#     'faxtor006a',
#     outputfilename=f'{path}/faxtor006a_highlights',
#     method=preprocessing.stack_movie.NONE,
#     output_type='frames',
#     figsize_mult = 1/3,
#     highlights=True
# )

  

for file in [
    # 'faxtor006a',
    # 'faxtor006a_hpf',
    # 'faxtor006a_hpf_movavrem',
    # 'faxtor006a_movavrem'
    # 'faxtor030a_hpf_movavrem',
    # 'faxtor030b_hpf_movavrem',
    # 'faxtor034_hpf_movavrem',
    # 'faxtor035_hpf_movavrem',
]:

    preprocessing.stack_movie.go(
        file,
        outputfilename=f'{path}/{file}',
        method=preprocessing.stack_movie.NONE,
        output_type='frames',
        figsize_mult = 1/3,
    )

# for file, fit in [
#     ('faxtor006a',              True),
#     ('faxtor006a_hpf',          True),
#     ('faxtor006a_hpf_movavrem', False),
#     ('faxtor006a_movavrem',     False),
# ]:
    
#     fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.4))
#     DDM.show_single.go(
#         file,
#         ax,
#         # target_ks=[0.37],
#         target_ks=[1.9],
#         fit_use_flow=True,
#         do_fit=fit,
#     )
#     common.save_fig(fig, f'{path}/ddm_overlapped_{file}.png', dpi=200)


# anatomix faxtor comparison

preprocessing.stack_movie.go(
    'faxtor030b_hpf_movavrem',
    outputfilename=f'{path}/faxtor030a_crop',
    method=preprocessing.stack_movie.NONE,
    output_type='frames',
    # figsize_mult = 1/3,
    crop = ((380, 380+420), (1307,1307+420)),
)
preprocessing.stack_movie.go(
    'pierre_exp',
    outputfilename=f'{path}/pierre_exp',
    method=preprocessing.stack_movie.NONE,
    output_type='frames',
    # figsize_mult = 1/3,
)