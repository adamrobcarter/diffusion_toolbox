# this file is just for seeing what's saved inside a data fi

import common

if __name__ == '__main__':
    for file in common.files_from_argv('particle_linking/data/', 'trajs_'):
        data = common.load(f'particle_linking/data/trajs_{file}.npz')
        particles = data['particles']
        diameter = data['particle_diameter']
        print(f'mean num particles, {particles.shape[0]/(particles[:, 2].max()+1):.0f}')
        print(f'{particles[:, 0].min():.3f} <= x <= {particles[:, 0].max():.3f}')
        print(f'{particles[:, 1].min():.3f} <= y <= {particles[:, 1].max():.3f}')
        if data.get('dimension', 2) == 3:
            print(f'{particles[:, 2].min():.3f} <= z <= {particles[:, 2].max():.3f}')
            print(f'<z> = {particles[:, 2].mean()/(diameter/2):.3f}a')