import common
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.patches
import tqdm

import box_counting.example

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    

    data_particles = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data_particles['particles']

    data_stack = common.load(f'preprocessing/data/stack_{file}.npz')
    stack = data_stack['stack']
    pixel_size = data_stack['pixel_size']

    FRAME = 0

    box_size_um = 37
    sep_size_um = 11

    every_nth_frame = 25
    for i, frame in enumerate(tqdm.trange(0, every_nth_frame*50, every_nth_frame)):

        box_counting.example.show(box_size_um, sep_size_um, f'box_counting/figures_png/example_{file}/example_{file}_{i}.png', True, particles, stack, pixel_size, frame,
                                  dpi=150)
