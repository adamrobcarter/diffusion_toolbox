import common
import numpy as np

CROPS = [800]

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    for crop in CROPS:
        size_x = crop
        size_y = crop
        start_x = (data['window_size_x'] - size_x) / 2
        start_y = (data['window_size_y'] - size_y) / 2
        print(size_x/data['window_size_x'], size_y/data['window_size_y'], start_x/data['window_size_x'], start_y/data['window_size_y'])
        particles_new = common.crop_particles(particles, start_x+size_x, start_y+size_y, start_x, start_y)
        ratio = particles_new.size / particles.size
        print(f'kept {ratio:.3f}')
        # assert np.isclose(crop**2, ratio, atol=0.1)
        newdata = dict(data)
        newdata['particles'] = particles_new
        newdata['window_size_x'] = size_x
        newdata['window_size_y'] = size_y
        common.save_data(f'particle_detection/data/particles_{file}_crop{crop}.npz', **newdata)
        print()