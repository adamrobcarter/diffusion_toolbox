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

FILES = ['sim_nohydro_011_L640_pot', 'sim_nohydro_002_L640_pot']

# for file in FILES:
#     particle_linking.link.go(f'{file}_div64')
#     MSD.calc.go(f'{file}_div64')
#     MSD.show.go(f'{file}_div64')

# for file in FILES:
#     particle_detection.crop.go(f'{file}',        crop=0.9)
#     particle_detection.crop.go(f'{file}_longer', crop=0.9)

# for file in FILES:
#     for window in True, False:
#         isf.calc_f.go(f'{file}_crop0.9',        window=window)
#         isf.calc_f.go(f'{file}_longer_crop0.9', window=window)

for file in FILES:
    isf.show_f.go(f'{file}_crop0.9')
    isf.show_f.go(f'{file}_longer_crop0.9')
    isf.show_f.go(f'{file}_crop0.9_bhwindow')
    isf.show_f.go(f'{file}_longer_crop0.9_bhwindow')

# for file in FILES:
#     for source in ['f_']
#     visualisation.merge_Ds

for file in FILES:

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    visualisation.Ds_overlapped_mult.go(
        files = [f'{file}_crop0.9', f'{file}_longer_crop0.9', f'{file}_crop0.9_bhwindow', f'{file}_longer_crop0.9_bhwindow'],
        ax = ax,
        legend_fontsize = 8,
        # linestyles = ['-'],
        colors = ['tab:blue', 'tab:blue', 'tab:orange', 'tab:orange'],
        sources = ['f_t0.5, f_t256'],
    )
    common.save_fig(fig, f'workflows/figures/windowing_{file}.png')