import preprocessing.stack_movie
import DDM.show_single
import common
import matplotlib.pyplot as plt

path = '/home/acarter/presentations/soil_24_dec/figures'

# preprocessing.stack_movie.go(
#     'psiche052',
#     f'{path}/psiche052.mp4',
#     method=preprocessing.stack_movie.REMOVE_BACKGROUND,
# )
# preprocessing.stack_movie.go(
#     'psiche053',
#     f'{path}/psiche053.mp4',
#     method=preprocessing.stack_movie.REMOVE_BACKGROUND,
# )

# fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(8, 4))

# ddm_kwargs = dict(
#     fit_result_in_label = True,
#     target_ks           = (0.125, 0.5, 1, 2, 4),
# )
# DDM.show_single.go('psiche052', ax0, **ddm_kwargs)
# DDM.show_single.go('psiche053', ax1, **ddm_kwargs)
# XLIM = (0.5, 600)
# YLIM = (1.2e3, 10e3)
# ax0.set_xlim(*XLIM)
# ax1.set_xlim(*XLIM)
# ax0.set_ylim(*YLIM)
# ax1.set_ylim(*YLIM)
# ax0.semilogy()
# ax1.semilogy()

# common.save_fig(fig, f'{path}/ddm.png', dpi=200)




# fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(8, 4))

# ddm_kwargs = dict(
#     fit_result_in_label = True,
#     target_ks           = (0.125, 0.5, 1, 2, 4),
# )
# DDM.show_single.go('psiche103', ax0, **ddm_kwargs)
# DDM.show_single.go('psiche103_shifted', ax1, **ddm_kwargs)
# XLIM = (0.5, 600)
# YLIM = (1.2e3, 10e5)
# # ax0.set_xlim(*XLIM)
# # ax1.set_xlim(*XLIM)
# ax0.set_ylim(*YLIM)
# ax1.set_ylim(*YLIM)
# ax0.semilogy()
# ax1.semilogy()

# common.save_fig(fig, f'{path}/ddm2.png', dpi=200)



preprocessing.stack_movie.go(
    'psiche103',
    f'{path}/psiche103.mp4',
    # method=preprocessing.stack_movie.REMOVE_BACKGROUND,
)