import common

TRIMS = [7.11]
# TRIMS = [0.0625]

for file in common.files_from_argv('particle_detection/data', 'particles_'):
    data = common.load(f'particle_detection/data/particles_{file}.npz')
    particles = data['particles']
    time_step = data['time_step']


    for trim in TRIMS:
        max_time_seconds = trim * 60 * 60
        max_time_frames = int(max_time_seconds / time_step)

        particles_in_crop = particles[:, 2] < max_time_frames
        newdata = dict(data)
        newdata['particles'] = particles[particles_in_crop, :]
        newdata['max_time_hours'] = trim
        common.save_data(f'particle_detection/data/particles_{file}_trim{trim}hr.npz', **newdata)