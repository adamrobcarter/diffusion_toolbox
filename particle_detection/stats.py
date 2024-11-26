# this file is just for seeing what's saved inside a data fi

import common

for file in common.files_from_argv('particle_detection/data/', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    print(f'mean num particles, {particles.shape[0]/(particles[:, 2].max()+1):.0f}')