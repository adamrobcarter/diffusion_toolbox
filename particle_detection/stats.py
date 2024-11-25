# this file is just for seeing what's saved inside a data fi

import common

for file in common.files_from_argv('particle_detection/data/', 'particles_'):
    common.load(f'particle_detection/data/particles_{file}.npz')