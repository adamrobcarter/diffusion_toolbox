import isf.calc_f
import isf.show_f
import visualisation.merge_Ds
import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import particle_detection.crop
import common
import particle_linking.link
import MSD.calc
import MSD.show
import preprocessing.brennan

# preprocessing.brennan.go_mesu('/store/cartera/2d_monolayer/nohydro_t16_pot_phi0.114_L640.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)
# preprocessing.brennan.go_mesu('/store/cartera/2d_monolayer/nohydro_t16_longest_pot_phi0.016_L640.bin', suffix='_pot_longest', multiply_time_by=1/16, dt=16)

# FILES = ['sim_nohydro_011_L640_pot_longer', 'sim_nohydro_002_L640_pot_longest']
FILES = ['sim_nohydro_011_L640_pot', 'sim_nohydro_002_L640_pot']

# for file in FILES:
#     particle_linking.link.go(f'{file}_div64')
#     MSD.calc.go(f'{file}_div64')
#     MSD.show.go(f'{file}_div64')

# for file in FILES:
#     particle_detection.crop.go(f'{file}', crop=0.9)

# for file in FILES:
#     for window in True, False:
#         isf.calc_f.go(f'{file}_crop0.9', window=window)

# for file in FILES:
#     isf.show_f.go(f'{file}_crop0.9')
#     isf.show_f.go(f'{file}_crop0.9_bhwindow')

# # for file in FILES:
# #    for source in ['f_']
# #    visualisation.merge_Ds

for file in ['sim_nohydro_011_L640_pot_longer', 'sim_nohydro_002_L640_pot_longest']:

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    visualisation.Ds_overlapped_mult.go(
        [
            (f'{file}_crop0.9',          'D0Sk_theory'),
            # (f'{file}_crop0.9',                 'f_t0.5'),
            (f'{file}_crop0.9',          'f_first'),
            # (f'{file}_crop0.9',          'f_t256'),
            # (f'{file}_crop0.9_bhwindow',        'f_t0.5'),
            (f'{file}_crop0.9_bhwindow',          'f_first'),
            # (f'{file}_crop0.9_bhwindow',          'f_t256'),
        ],
        ax = ax,
        legend_fontsize = 8,
        # linestyles = ['-'],
        discrete_colors=True
    )
    common.save_fig(fig, f'workflows/figures/windowing_{file}.png')

for file in ['sim_nohydro_011_L640_pot', 'sim_nohydro_002_L640_pot']:

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    visualisation.Ds_overlapped_mult.go(
        [
            (f'{file}_crop0.9',          'D0Sk_theory'),
            (f'{file}_crop0.9',                 'f_t0.5'),
            # (f'{file}_crop0.9',          'f_first'),
            (f'{file}_crop0.9',          'f_t128'),
        ],
        ax = ax,
        legend_fontsize = 8,
        # linestyles = ['-'],
        discrete_colors=True,
    )
    ax.set_ylim(0.7, 6)
    common.save_fig(fig, f'workflows/figures/windowing_fkt_divergence_time_{file}.png')