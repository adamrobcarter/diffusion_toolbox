# this file is just for seeing what's saved inside a data fi

import common

for file in common.files_from_argv('particle_detection/data/', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    print(f'mean num particles, {particles.shape[0]/(particles[:, 2].max()+1):.0f}')
    print(f'{particles[:, 0].min()} <= x <= {particles[:, 0].max()}')
    print(f'{particles[:, 1].min()} <= y <= {particles[:, 1].max()}')