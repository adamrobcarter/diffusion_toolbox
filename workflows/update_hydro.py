import isf.calc_f
import box_counting.count
import isf.show_f
import box_counting.msd_single
import box_counting.D_of_L
import visualisation.merge_Ds
import visualisation.Ds_overlapped_mult
import matplotlib.pyplot as plt
import common
import preprocessing.brennan

# preprocessing.brennan.go_mesu('/store/cartera/2d_monolayer/hydro_t16_pot_phi0.114_L320.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)
preprocessing.brennan.go_mesu('/store/cartera/2d_monolayer/hydro_t16_pot_phi0.114_L640.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)

FILES = ['sim_hydro_011_L320_pot_longer', 'sim_hydro_011_L640_pot_longer']
FILES = ['sim_hydro_011_L640_pot_longer']

for file in FILES:
    isf.calc_f.go(file)

for file in FILES:
    box_counting.count.go(file)

for file in FILES:
    isf.show_f.go(file, skip_saving_figure=True)

    
for file in FILES:
    box_counting.msd_single.go(file)
    
for file in FILES:
    box_counting.D_of_L.go(file, plateau_source='var')
    
SOURCES = ['timescaleint_nofit_cropped_var', 'timescaleint_fixexponent_var', 'f_t64', 'f_t256', 'f_t1024']
# for source in SOURCES:
#     visualisation.merge_Ds.go(['sim_hydro_011_L320_pot', 'sim_hydro_011_L320_pot_longer'], source)
#     visualisation.merge_Ds.go(['sim_hydro_011_L640_pot', 'sim_hydro_011_L640_pot_longer'], source)

D_of_L_props = dict(
    legend_fontsize = 8,
    linestyles = ['-']*len(SOURCES),
    colors = ['tab:blue', 'tab:cyan', common.colormap(0.25), common.colormap(0.5), common.colormap(0.75)],
)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
visualisation.Ds_overlapped_mult.go(
    [('sim_hydro_011_L320_pot_longer', source) for source in SOURCES],
    ax = ax,
    **D_of_L_props
)
common.save_fig(fig, 'workflows/hydro_L320.png')

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
visualisation.Ds_overlapped_mult.go(
    [('sim_hydro_011_L640_pot_longer', source) for source in SOURCES],
    ax = ax,
    **D_of_L_props
)
common.save_fig(fig, 'workflows/hydro_L640.png')