import common

TRIMS = [1.0, 0.5, 0.25, 0.125, 0.0625]
# TRIMS = [0.0625]

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    for trim in TRIMS:
        particles_in_crop = particles[:, 2] < particles[:, 2].max() * trim
        newdata = dict(data)
        newdata['particles'] = particles[particles_in_crop, :]
        newdata['max_time_hours'] = data['max_time_hours'] * trim
        common.save_data(f'particle_detection/data/particles_{file}_trim{trim}.npz', **newdata)