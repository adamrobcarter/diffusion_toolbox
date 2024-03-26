import numpy as np
import sys

import box_counting.count

for file in sys.argv[1:]:
    box_sizes_px = np.logspace(1, 9, 20, base=2)
    sep_sizes_px = np.array([25]*box_sizes_px.shape[0])

    output_file_name = f'box_counting/data/counted_dense_{file}.npz'

    box_counting.count.calc_and_save(box_sizes_px, sep_sizes_px, file, output_file_name)