import numpy as np
import sys
import matplotlib.pyplot as plt

file = sys.argv[1]

data = np.load(f'intensity_counting/data/all_counts_{file}.npz')
counts_intensity = data['counts']
box_sizes_intensity = data['box_sizes']

data = np.load(f'box_counting/data/all_counts_{file}.npz')
counts_boxes = data['counts']
box_sizes_boxes = data['box_sizes']
print(counts_boxes.shape, counts_intensity.shape)

fig, axs = plt.subplots(10, 4, figsize=(14, 4*10))

if file == 'pierre_sim':
    intensity_factor = 16.17**2
counts_intensity /= np.sqrt(intensity_factor)

boxes_to_use = (0, 2, 4, 6)
ymaxs = np.zeros(len(boxes_to_use))
ymins = np.zeros(len(boxes_to_use))

for chosen_box_index in range(10):
    for box_size_index_index, box_size_index in enumerate(boxes_to_use):
        ax = axs[chosen_box_index][box_size_index_index]
        
        N_of_t = counts_boxes    [box_size_index, chosen_box_index, :]
        N_of_0 = counts_boxes    [box_size_index, chosen_box_index, 0]
        delta_N = N_of_t - N_of_0
        I_of_t = counts_intensity[box_size_index, :, chosen_box_index]
        I_of_0 = counts_intensity[box_size_index, 0, chosen_box_index]
        delta_I = I_of_t - I_of_0
        ax.plot(delta_N, label=f'boxes {    N_of_t.var():.3f}', zorder=3, alpha=0.7)
        ax.plot(delta_I, label=f'intensity {I_of_t.var():.3f}')
        ax.set_title(f'{box_sizes_boxes[box_size_index]} {box_sizes_intensity[box_size_index]}')
        ax.legend()

        max_ = max(np.max(delta_N), np.max(delta_I))
        min_ = min(np.min(delta_N), np.min(delta_I))
        ymaxs[box_size_index_index] = max(max_, ymaxs[box_size_index_index])
        ymins[box_size_index_index] = min(min_, ymins[box_size_index_index])

# set all ylims to the same
for box_size_index_index in range(len(boxes_to_use)):
    border = 0.03*(ymaxs[box_size_index_index] - ymins[box_size_index_index])
    for ax in axs[:, box_size_index_index]:
        ax.set_ylim(ymins[box_size_index_index]-border, ymaxs[box_size_index_index]+border)

fig.tight_layout()
fig.savefig(f'visualisation/figures_png/individual_box_counts_{file}.png', dpi=200)