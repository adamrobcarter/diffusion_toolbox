import common
import numpy as np

# CROPS = [1.0, 0.5, 0.25, 0.125, 0.0625]
CROPS = [0.9]

def go(file, crop):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    print(f'{particles[:, 0].min()} <= x <= {particles[:, 0].max()}')
    print(f'{particles[:, 1].min()} <= y <= {particles[:, 1].max()}')

    for crop in CROPS:
        size_x = data['window_size_x'] * crop
        size_y = data['window_size_y'] * crop
        if '_pot' in file and not 'nohydro' in file: # it seems only hydro pot goes [-L/2, L/2], nohydro pot goes [0, L]
            start_x = - data['window_size_x'] * crop / 2
            start_y = - data['window_size_y'] * crop / 2
        else:
            start_x = (data['window_size_x'] - size_x) / 2
            start_y = (data['window_size_y'] - size_y) / 2
        print(size_x/data['window_size_x'], size_y/data['window_size_y'], start_x/data['window_size_x'], start_y/data['window_size_y'])
        particles_new = common.crop_particles(particles, start_x+size_x, start_y+size_y, start_x, start_y)
        ratio = particles_new.size / particles.size
        print(f'kept {ratio:.3f}')
        assert np.isclose(crop**2, ratio, atol=0.1)
        newdata = dict(data)
        newdata['particles'] = particles_new
        newdata['window_size_x'] = size_x
        newdata['window_size_y'] = size_y
        common.save_data(f'particle_detection/data/particles_{file}_crop{crop}.npz', **newdata)
        print()

if __name__ == '__main__':
    for file in common.files_from_argv('particle_detection/data', 'particles_'):
        for crop in CROPS:
            go(file, crop)