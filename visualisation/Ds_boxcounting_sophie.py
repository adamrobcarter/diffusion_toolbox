import common
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)

timescaleintegral_x = np.array([	0.103942652329749,0.207885304659498,0.519713261648746,1.03942652329749,1.43440860215054,2.07885304659498,2.87921146953405,4.15770609318996,5.75842293906810,8.31541218637993,11.5064516129032,16.6308243727599,19.7491039426523,22.8673835125448,24.4265232974910])
timescaleintegral_y = np.array([0.715084639456640,0.759462918744972,1.46261640042440,2.56073086204352,3.78909538953524,5.41105358248309,6.78255486238759,8.58502158959043,8.95193674486290,9.93854195925673,9.47657037534071,9.06689661448547,9.39221079633647,9.01776477004307,8.02304202449399])

data_msd = common.load(f'visualisation/data/Ds_from_MSD_short_eleanor0.01')
D_msd = data_msd['Ds'][0]

data_tsi = common.load(f'visualisation/data/Ds_from_timescaleint_eleanorlong')
Ds_tsi = data_tsi['Ds']
sigma  = data_tsi['particle_diameter']
L_tsi  = data_tsi['Ls']/sigma
ax.scatter(L_tsi, Ds_tsi/D_msd, label='Adam timescale integral')

data_tsi_nofit = common.load(f'visualisation/data/Ds_from_timescaleint_nofit_eleanorlong')
Ds_tsi_nofit = data_tsi_nofit['Ds']
L_tsi_nofit  = data_tsi_nofit['Ls']/sigma
ax.scatter(L_tsi_nofit, Ds_tsi_nofit/D_msd, label='Adam timescale integral nofit')

data_fit = common.load(f'visualisation/data/Ds_from_boxcounting_collective_eleanorlong')
Ds_fit = data_fit['Ds']
L_fit  = data_fit['Ls']/sigma
ax.scatter(L_fit, Ds_fit/D_msd, label='Adam NMSD fit')

ax.scatter(timescaleintegral_x, timescaleintegral_y, label='from Sophie', marker='x')

ax.legend(fontsize=7)
ax.loglog()
ax.set_ylim(0.5, 70)
common.save_fig(fig, 'visualisation/figures_png/Ds_boxcounting_sophie.png')

