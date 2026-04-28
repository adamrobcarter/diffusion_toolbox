import van_hove.show_g_mult
import matplotlib.pyplot as plt
import workflows.thesis.common
import common

files = ('eleanorlong034', 'eleanorlong010')

# fig, ax = plt.subplots(figsize=workflows.thesis.common.figsize_small)

# for i, file in enumerate(files):
#     # color = plt.cm.copper(i/(len(files))-0.5)
#     van_hove.show_g_mult.go(file, ax)

# ax.set_xlim(0, 10)

# ax.legend()
# ax.axhline(1.0, color='gray', linestyle=':')
# common.save_fig(fig, f'workflows/thesis/figures/g_of_r.pdf', hide_metadata=True)


import isf.show_S_of_k

fig, ax = plt.subplots(figsize=workflows.thesis.common.figsize_small)

for i, file in enumerate(files):
    print('len(files)', len(files))
    color = plt.cm.copper(i/(len(files)-0.5))
    print(i, i/(len(files)-0.5), 'color', color)
    isf.show_S_of_k.go(file, ax, source='F_first', theory_color='gray', data_color=color, show_fit=True)

# ax.set_xlim(0, 10)

ax.legend()
ax.axhline(1.0, color='gray', linestyle=':')
common.save_fig(fig, f'workflows/thesis/figures/S_of_k.pdf', hide_metadata=True)