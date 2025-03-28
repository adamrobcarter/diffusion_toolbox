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
# preprocessing.brennan.go_mesu('/store/cartera/2d_monolayer/hydro_t16_pot_phi0.114_L640.bin', suffix='_pot_longer', multiply_time_by=1/16, dt=16)

isf.calc_f.go('sim_hydro_011_L320_pot_longer')
isf.calc_f.go('sim_hydro_011_L640_pot_longer')
box_counting.count.go('sim_hydro_011_L320_pot_longer')
box_counting.count.go('sim_hydro_011_L640_pot_longer')

isf.show_f.go('sim_hydro_011_L320_pot_longer', skip_saving_figure=True)
isf.show_f.go('sim_hydro_011_L640_pot_longer', skip_saving_figure=True)
box_counting.msd_single.go('sim_hydro_011_L320_pot_longer')
box_counting.msd_single.go('sim_hydro_011_L640_pot_longer')
box_counting.D_of_L.go('sim_hydro_011_L320_pot_longer', plateau_source='var')
box_counting.D_of_L.go('sim_hydro_011_L640_pot_longer', plateau_source='var')
SOURCES = ['timescaleint_nofit_cropped_var', 'f_t64', 'f_t256', 'f_t1024']
# for source in SOURCES:
#     visualisation.merge_Ds.go(['sim_hydro_011_L320_pot', 'sim_hydro_011_L320_pot_longer'], source)
#     visualisation.merge_Ds.go(['sim_hydro_011_L640_pot', 'sim_hydro_011_L640_pot_longer'], source)

D_of_L_props = dict(
    legend_fontsize = 8,
    linestyles = ['-'],
    colors = [['tab:blue', common.colormap(0.25), common.colormap(0.5), common.colormap(0.75)]],
    sources = SOURCES,
)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
visualisation.Ds_overlapped_mult.go(
    files = ['sim_hydro_011_L320_pot_longer'],
    ax = ax,
    **D_of_L_props
)
common.save_fig(fig, 'workflows/hydro_L320.png')

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
visualisation.Ds_overlapped_mult.go(
    files = ['sim_hydro_011_L640_pot_longer'],
    ax = ax,
    **D_of_L_props
)
common.save_fig(fig, 'workflows/hydro_L640.png')