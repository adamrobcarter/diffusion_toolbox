import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches

box_size_x = 128 * 0.17
box_size_y = 128 * 0.17
sep_size = -40 * 0.17
# sep_size = 20 * 0.17
window_size_x=214.3
window_size_y=170.8

SepSize_x = box_size_x + sep_size
SepSize_y = box_size_y + sep_size
num_boxes_x = int(np.floor(window_size_x / SepSize_x))
num_boxes_y = int(np.floor(window_size_y / SepSize_y))
print(num_boxes_x, num_boxes_y)

fig, ax = plt.subplots(1, 1)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

i = 0

for box_x_index in range(num_boxes_x):
    for box_y_index in range(num_boxes_y):

        box_x_min = box_x_index * SepSize_x# + sep_sizes[box_index]/2
        box_x_max = box_x_min + box_size_x
        box_y_min = box_y_index * SepSize_y# + sep_sizes[box_index]/2
        box_y_max = box_y_min + box_size_y

        if i > len(colors)-1:
            i = 0
        color = colors[i]
        i += 1

        rect = matplotlib.patches.Rectangle((box_x_min, box_y_min), box_size_x, box_size_y,
                                            linewidth=1, edgecolor='black', facecolor='none')#, alpha=0.05)

        ax.add_patch(rect)

# rect = matplotlib.patches.Rectangle((50,50), 50, 30, linewidth=1, edgecolor='b', facecolor='r')
# ax.add_patch(rect)
# ax.set_xlim(left = 0, right = 150)
# ax.set_ylim(bottom = 0, top = 150)

ax.set_xlim(0, window_size_x)
ax.set_ylim(0, window_size_y)
fig.savefig('box_counting/figures_png/shown_boxes.png', dpi=300)