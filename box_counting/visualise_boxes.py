import common
import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np
import matplotlib.cm

for file in common.files_from_argv('box_counting/data', 'counted_'):
    data = common.load(f'box_counting/data/counted_{file}.npz')
    

    box_sizes = data['box_sizes']
    # box_size_index = len(box_sizes)//2 + 1
    box_size_index = 0
    box_coords = data['box_coords'][box_size_index, :, :, :]

    xs = box_coords[:, :, 0]
    ys = box_coords[:, :, 1]

    fig, ax = plt.subplots(1, 1)

    for i in range(xs.shape[0]):
        for j in range(xs.shape[1]):
            color = matplotlib.cm.hsv(((i*xs.shape[0] + j)*np.pi*10)%1)
            rect = matplotlib.patches.Rectangle((xs[i, j], ys[i, j]), box_sizes[box_size_index], box_sizes[box_size_index],
                linewidth=1, facecolor=color, alpha=0.4)
            ax.add_patch(rect)

    ax.set_aspect('equal')
    border = data['window_size_x'] / 10
    ax.set_xlim(-border, data['window_size_x']+border)
    ax.set_ylim(-border, data['window_size_y']+border)

    rect = matplotlib.patches.Rectangle((0, 0), data['window_size_x'], data['window_size_y'],
        linewidth=1, edgecolor='black', facecolor='none')
    ax.add_patch(rect)

    common.save_fig(fig, f'box_counting/figures_png/visualise_boxes_{file}.png', dpi=300)