import common

CROP = 0.25

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    crop_x = data['window_size_x'] * CROP
    crop_y = data['window_size_y'] * CROP
    particles_in_crop = (particles[:, 0] < crop_x) & (particles[:, 1] < crop_y)
    newdata = dict(data)
    newdata['particles'] = particles[particles_in_crop, :]
    newdata['window_size_x'] = crop_x
    newdata['window_size_y'] = crop_y
    common.save_data(f'particle_detection/data/particles_{file}_crop{CROP}.npz', **newdata)