import common
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)

for file in (files := common.files_from_argv('visualisation/data', 'Ds_from_f')):
    data = common.load(f'visualisation/data/Ds_from_f_{file}.npz')
    ks     = data['ks']
    Ds     = data['Ds']
    D_uncs = data['D_uncs']
    particle_diameter = data['particle_diameter']
    
    data_msd = common.load(f'visualisation/data/Ds_from_MSD_{file}.npz')
    D_MSD     = data_msd['Ds'][0]
    pack_frac = data_msd['pack_frac_given']

    lambdas = 2 * np.pi / ks

    Lh = 2 / ( 3 * pack_frac ) * particle_diameter
    Lh = 1

    ax.errorbar(Lh/lambdas, Ds/D_MSD, yerr=D_uncs/D_MSD, label=file, linestyle='none', marker='o')

x_th = np.logspace(-1, 1)
ax.plot(x_th, 1+1/(2*np.pi)/x_th, color='gray')

ax.set_xlabel(r'$2\pi/k / L_h$')
ax.set_ylabel(r'$D(k)/D_{MSD}$')
ax.loglog()

filenames = '_'.join(files)
common.save_fig(fig, f'visualisation/figures_png/panzuela_recaling_{filenames}.png')