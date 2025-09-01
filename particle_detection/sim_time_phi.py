import common
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

phis = [0.02, 0.04, 0.06, 0.08, 0.1]
filegroups = [
    [f'sim_hydro_{phi}_L500_t1e2_0.5' for phi in phis],
    [f'sim_hydro_{phi}_L500_t1e2_0.5_pot' for phi in phis],
]
filegroup_names = ['DPStokes (p, p, sw)', 'NBody (o, o, sw)']

for filegroup_i, filegroup in enumerate(filegroups):
    phis = []
    ts = []

    print('filegroup', filegroup)

    for file in filegroup:
        # print('files', files)
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        computation_time = data['computation_time']
        dt = data.get('dt', 0.1)
        saved_time_step = data['time_step']
        time_column = data['dimension']
        L  =data['window_size_x']
        particle_diameter = data['particle_diameter']
        num_saved_timesteps = np.unique(data['particles'][:, time_column]).size
        max_time = num_saved_timesteps * saved_time_step
        num_sim_timesteps = max_time / dt
        print('num timesteps', num_sim_timesteps)
        time_per_step = computation_time / num_sim_timesteps
        # phis.append(data['pack_frac'])
        N = data['density'] * data['window_size_x'] * data['window_size_y']
        N_over_phi = N/data['pack_frac']
        ts.append(time_per_step*1000)
        phis.append(N)

    popt, unc = common.curve_fit(lambda x, m, c: m*x + c, phis, ts)
    label = filegroup_names[filegroup_i]
    label += f' $t={popt[0]:.4f} \phi + {popt[1]:.0f}$'
    ax.plot(phis, ts, marker='o', linestyle=':', label=label)
    

ax.text(0.05, 0.8, f'$L={L/particle_diameter:.0f}\sigma$', transform=ax.transAxes, fontsize=12)

# ax.set_xlabel(f'packing fraction $\phi$ ($N={N_over_phi:.0f}\phi$)')
ax.set_xlabel('number of particles')
ax.set_ylabel('time per simulation step (ms)')
ax.legend(fontsize=10)
# ax.grid()

# common.add_exponential_index_indicator(ax, 2, (6e2, 1e-2), 'L')

filestr = '_'.join(['_'.join(files) for files in filegroups])[:50]
common.save_fig(fig, f'particle_detection/figures_png/sim_time_phi_{filestr}.png')