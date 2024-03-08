import matplotlib.pyplot as plt
import numpy as np
import common
import sys

file = sys.argv[1]

for file in sys.argv[1:]:
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

    data = common.load(f'intensity_counting/data/counted_{file}.npz')
    box_sizes               = data['box_sizes']
    counted_intensity_diffs = data['counted_intensity_diffs']
    avg_intensities         = data['avg_intensities']
    pixel_size              = data['pixel_size']
    time_step               = data['time_step']
    variances               = data['variances']
    particle_diameter       = data['particle_diameter']

    plt.scatter(box_sizes, variances/box_sizes**2)
    plt.loglog()
    plt.xlabel('L')
    plt.ylabel('plateau / L^2')
    plt.tight_layout()
    plt.savefig(f'visualisation/figures_png/var_vs_L_{file}.png', dpi=300)
