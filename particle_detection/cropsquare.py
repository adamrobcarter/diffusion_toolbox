import common
import numpy as np

CROP = np.sqrt(2)/2
CROP = 1
CROP = 0.9

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    crop_width = min(data['window_size_x'], data['window_size_y']) * CROP
    crop_start = (min(data['window_size_x'], data['window_size_y']) - crop_width) / 2
    particles_new = common.crop_particles(particles, crop_start+crop_width, crop_start+crop_width, crop_start, crop_start)
    newdata = dict(data)
    newdata['particles'] = particles_new
    newdata['window_size_x'] = crop_width
    newdata['window_size_y'] = crop_width
    common.save_data(f'particle_detection/data/particles_{file}_cropsquare{CROP}.npz', **newdata)