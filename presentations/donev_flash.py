import common
import matplotlib.pyplot as plt

path = '/home/acarter/presentations/donev_flash/figures'

import particle_detection.show_movie
particle_detection.show_movie.go(
    'ld_hydro_nbody_0.07_singlewall_L160_t1h_1_uneven',
    # infile = f'particle_detection/data/particles_{file}.npz',
    outfile = f'/home/acarter/presentations/donev_flash/figures/tophat',
    output_type='frames',
    dpi=100,
)

import box_counting.show_equilib
box_counting.show_equilib.go(
    file = 'ld_hydro_nbody_0.07_singlewall_L160_t75s_1',
    force_fps = 10,
    outputfilename = '/home/acarter/presentations/donev_flash/figures/equilib',
    output_type = 'frames',
    dpi=100,
)

# import visualisation.Ds_overlapped_mult

# fig, ax = plt.subplots(figsize=(4, 2.8))

# visualisation.Ds_overlapped_mult.go(
#     [
#         dict(
#             file = 'sim_nohydro_011_L1280_longer_mergedD',
#             source = 'f_first_first',
#             marker = 'o',
#         )
#     ],
#     ax = ax,
#     plot_against_k = True,
#     show_twin_k_axis = False,
#     allow_rescale_x = True,
# )
# ax.get_legend().remove()
# ax.set_xlim(15e-2, 4e1)
# ax.set_ylim(0.85, 1.7)

# common.save_fig(fig, f'{path}/Dk_nohydro.pdf', hide_metadata=True)



import particle_detection.show_movie
particle_detection.show_movie.go(
    'ld_hydro_nbody_0.01_singlewall_L160_t75s_1',
    infile = f'particle_linking/data/trajs_ld_hydro_nbody_0.01_singlewall_L160_t10m_8s.npz',
    outfile = f'{path}/linked',
    # crop = 100,
    # every_nth_frame=10
    tracks = True,
    highlights = True,
    output_type='frames',
    show_blobs_too=True,
    dpi=100,
)