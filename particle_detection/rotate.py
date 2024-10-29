import common

ROTATION = 45

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    width     = data['window_size_x']
    height    = data['window_size_y']

    particles_new, width_new, height_new = common.rotate_particles(ROTATION, particles, width, height)
    
    newdata = dict(data)
    newdata['particles']     = particles_new
    newdata['window_size_x'] = width_new
    newdata['window_size_y'] = height_new
    common.save_data(f'particle_detection/data/particles_{file}_rot{ROTATION}.npz', **newdata)