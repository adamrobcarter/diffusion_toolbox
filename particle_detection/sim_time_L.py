import common
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

# filegroups = [
#     [f'sim_hydro_0.114_L{L}_t1e3_0.5' for L in [160, 320, 640, 1280]],
#     [f'sim_hydro_0.114_L{L}_t1e3_0.5_pot' for L in [160, 320, 640, 1280]],
# ]
Ls = [int(L) for L in np.logspace(np.log10(100), np.log10(2000), num=10)]
# Ls = Ls
filegroups = [
    [f'sim_hydro_0.114_L{L}_t1e3_0.5' for L in Ls],
    [f'sim_hydro_0.114_L{L}_t1e3_0.5_pot' for L in Ls],
]
filegroup_names = ['DPStokes (p, p, sw)', 'NBody (o, o, sw)']

for filegroup_i, filegroup in enumerate(filegroups):
    Ls = []
    ts = []

    print('filegroup', filegroup)

    for file in filegroup:
        # print('files', files)
        data = common.load(f'particle_detection/data/particles_{file}.npz')
        computation_time = data['computation_time']
        dt = data.get('dt', 0.1)
        saved_time_step = data['time_step']
        time_column = data['dimension']
        pack_frac = data['pack_frac']
        num_saved_timesteps = np.unique(data['particles'][:, time_column]).size
        max_time = num_saved_timesteps * saved_time_step
        num_sim_timesteps = max_time / dt
        print('num timesteps', num_sim_timesteps)
        time_per_step = computation_time / num_sim_timesteps
        N = data['density'] * data['window_size_x'] * data['window_size_y']
        # Ls.append(data['window_size_x'])
        Ls.append(N)
        ts.append(time_per_step)

    ax.plot(Ls, ts, marker='o', linestyle=':', label=filegroup_names[filegroup_i])

ax.text(0.05, 0.5, f'$\phi={pack_frac:.3f}$', transform=ax.transAxes, fontsize=12)

ax.semilogx()
ax.semilogy()
ax.set_xlabel('simulation box size $L$')
ax.set_ylabel('time per simulation step (s)')
ax.legend()
# ax.grid()

common.add_exponential_index_indicator(ax, 2, (6e2, 1e-2), 'L')

filestr = '_'.join(['_'.join(files) for files in filegroups])[:50]
common.save_fig(fig, f'particle_detection/figures_png/sim_time_L_{filestr}.png')