import common

CROPS = [1.0, 0.5, 0.25, 0.125, 0.0625]

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    for crop in CROPS:
        crop_x = data['window_size_x'] * crop
        crop_y = data['window_size_y'] * crop
        particles_new = common.crop_particles(particles, crop_x, crop_y)
        newdata = dict(data)
        newdata['particles'] = particles_new
        newdata['window_size_x'] = crop_x
        newdata['window_size_y'] = crop_y
        common.save_data(f'particle_detection/data/particles_{file}_crop{crop}.npz', **newdata)