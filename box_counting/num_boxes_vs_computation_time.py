import common
import matplotlib.pyplot as plt
import sys
import numpy as np

NUM_SPLITS = 10
SPACINGS = [512.5, 256.5, 128.5, 64.5, 32.5, 16.5, 8.5, 4.5, 2.5, 1.5, 0.5]

for file in sys.argv[1:]:
    computation_times = []
    computation_time_stds = []
    num_of_boxes      = []

    for ax_i, spacing in enumerate(SPACINGS):

        individual_computation_times = []

        for split in range(NUM_SPLITS):
            data = common.load(f'box_counting/data/counted_{file}_qo_split{split}_spacing{spacing}.npz')
            
            N_stats = data['N_stats']
            individual_num_of_boxes = N_stats[:, 5].mean() # they should all be the same

            if individual_num_of_boxes < 100 and data['computation_time'] > 0.1:
                continue
            individual_computation_times.append(data['computation_time'])
            # print(individual_computation_times)
            # print(type(data['computation_time']))
            
        print(individual_computation_times)
        computation_times.append(np.mean(individual_computation_times))
        computation_time_stds.append(np.std(individual_computation_times))
        print(individual_num_of_boxes)
        num_of_boxes.append(individual_num_of_boxes)
        print()
        # breaks
    print(num_of_boxes)
    print(computation_times)
    plt.errorbar(num_of_boxes, computation_times, yerr=computation_time_stds, linestyle='none', marker='o')
    plt.semilogx()
    plt.semilogy()
    plt.xlabel('number of boxes')
    plt.ylabel('computation time (s)')
    common.save_fig(plt.gcf(), f'box_counting/figures_png/num_boxes_vs_computation_time_{file}.png')