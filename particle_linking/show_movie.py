

import numpy as np
import common
import sys
import particle_detection.show_movie


for file in sys.argv[1:]:
    particle_detection.show_movie.go(
        file,
        infile = f'particle_linking/data/trajs_{file}.npz',
        outfile = f'particle_linking/figures_png/movie_{file}.gif'
    )